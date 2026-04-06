from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from autodias.adapters import get_adapter
from autodias.config import load_config
from autodias.geometry import HARTREE_TO_KCAL, angle, dihedral, distance
from autodias.models import AnalysisRow, ProjectConfig, StructureSet
from autodias.structures import load_structures


def prepare(config_path: str | Path) -> None:
    config = load_config(config_path)
    structures = load_structures(config)
    _initialize_workspace(config)
    _write_log_header(config, structures)
    _write_xyz_files(config, structures)
    _write_input_files(config, structures)


def run(config_path: str | Path) -> None:
    config = load_config(config_path)
    structures = load_structures(config)
    _initialize_workspace(config)
    _write_log_header(config, structures)
    _write_xyz_files(config, structures)
    _write_input_files(config, structures)
    if config.prepare_only:
        _append_log(config, "No calculations requested.\n")
        return

    adapter = get_adapter(config.engine_name)
    for index in range(len(structures.frames)):
        _run_job(config, "complex", index)
        _run_job(config, config.fragment1_name, index)
        _run_job(config, config.fragment2_name, index)
        if config.analysis_enabled:
            row = _analyze_frame(config, structures, adapter, index)
            _append_analysis_row(config, row, initialize=(index == 0))
    _cleanup(config)


def analyze(config_path: str | Path) -> None:
    config = load_config(config_path)
    structures = load_structures(config)
    adapter = get_adapter(config.engine_name)
    for index in range(len(structures.frames)):
        row = _analyze_frame(config, structures, adapter, index)
        _append_analysis_row(config, row, initialize=(index == 0))


def _initialize_workspace(config: ProjectConfig) -> None:
    config.input_dir.mkdir(exist_ok=True)
    config.output_dir.mkdir(exist_ok=True)
    if config.keep_xyz:
        config.xyz_dir.mkdir(exist_ok=True)


def _write_log_header(config: ProjectConfig, structures: StructureSet) -> None:
    lines = [
        "autoDIAS 2.0\n",
        "-----------\n",
        f"Project: {config.name}\n",
        f"Input: {config.input_path}\n",
        f"Engine: {config.engine_name}\n",
        f"Frames: {len(structures.frames)}\n",
        f"Fragment 1 ({config.fragment1_name}): {structures.fragment1_atoms}\n",
        f"Fragment 2 ({config.fragment2_name}): {structures.fragment2_atoms}\n",
        "\n",
    ]
    config.log_path.write_text("".join(lines), encoding="utf-8")


def _append_log(config: ProjectConfig, message: str) -> None:
    with config.log_path.open("a", encoding="utf-8") as handle:
        handle.write(message)


def _write_xyz_files(config: ProjectConfig, structures: StructureSet) -> None:
    if not config.keep_xyz:
        return
    for index, (title, frame) in enumerate(zip(structures.titles, structures.frames, strict=True), start=1):
        _write_xyz(config.xyz_dir / f"complex_{index:04d}.xyz", structures.atoms, frame, title)
        _write_xyz(
            config.xyz_dir / f"{config.fragment1_name}_{index:04d}.xyz",
            _fragment_atoms(structures.atoms, structures.fragment1_atoms),
            _fragment_coords(frame, structures.fragment1_atoms),
            f"{config.fragment1_name}_{index:04d}",
        )
        _write_xyz(
            config.xyz_dir / f"{config.fragment2_name}_{index:04d}.xyz",
            _fragment_atoms(structures.atoms, structures.fragment2_atoms),
            _fragment_coords(frame, structures.fragment2_atoms),
            f"{config.fragment2_name}_{index:04d}",
        )


def _write_input_files(config: ProjectConfig, structures: StructureSet) -> None:
    adapter = get_adapter(config.engine_name)
    for index, frame in enumerate(structures.frames, start=1):
        complex_coords = _format_coordinates(structures.atoms, frame)
        frag1_coords = _format_coordinates(
            _fragment_atoms(structures.atoms, structures.fragment1_atoms),
            _fragment_coords(frame, structures.fragment1_atoms),
        )
        frag2_coords = _format_coordinates(
            _fragment_atoms(structures.atoms, structures.fragment2_atoms),
            _fragment_coords(frame, structures.fragment2_atoms),
        )

        (config.input_dir / f"complex_{index:04d}.{config.input_extension}").write_text(
            adapter.render_input(
                config.input_template,
                complex_coords,
                config.charge,
                config.multiplicity,
                f"complex_{index:04d}",
            ),
            encoding="utf-8",
        )
        (config.input_dir / f"{config.fragment1_name}_{index:04d}.{config.input_extension}").write_text(
            adapter.render_input(
                config.input_template,
                frag1_coords,
                config.fragment1_charge,
                config.fragment1_multiplicity,
                f"{config.fragment1_name}_{index:04d}",
            ),
            encoding="utf-8",
        )
        (config.input_dir / f"{config.fragment2_name}_{index:04d}.{config.input_extension}").write_text(
            adapter.render_input(
                config.input_template,
                frag2_coords,
                config.fragment2_charge,
                config.fragment2_multiplicity,
                f"{config.fragment2_name}_{index:04d}",
            ),
            encoding="utf-8",
        )


def _write_xyz(path: Path, atoms: list[str], coords: object, title: str) -> None:
    rows = [f"{len(atoms)}", title]
    for atom, xyz in zip(atoms, coords, strict=True):
        rows.append(f"{atom:2s} {xyz[0]:15.12f} {xyz[1]:15.12f} {xyz[2]:15.12f}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _format_coordinates(atoms: list[str], coords: object) -> str:
    return "\n".join(
        f"{atom:2s} {xyz[0]:15.12f} {xyz[1]:15.12f} {xyz[2]:15.12f}"
        for atom, xyz in zip(atoms, coords, strict=True)
    )


def _fragment_coords(frame, atom_indices: list[int]):
    return [frame[index - 1] for index in atom_indices]


def _fragment_atoms(atoms: list[str], atom_indices: list[int]) -> list[str]:
    return [atoms[index - 1] for index in atom_indices]


def _run_job(config: ProjectConfig, stem: str, index: int) -> None:
    input_path = config.input_dir / f"{stem}_{index + 1:04d}.{config.input_extension}"
    output_path = config.output_dir / f"{stem}_{index + 1:04d}.{config.output_extension}"
    command = config.command_template.replace("$input", str(input_path)).replace("$output", str(output_path))
    _append_log(config, f"Running {stem}_{index + 1:04d}: {command}\n")
    subprocess.run(command, shell=True, check=True, cwd=config.base_dir)


def _analyze_frame(config: ProjectConfig, structures: StructureSet, adapter, index: int) -> AnalysisRow:
    frame = structures.frames[index]
    complex_scf = adapter.parse_energy(config.output_dir / f"complex_{index + 1:04d}.{config.output_extension}")
    fragment1_scf = adapter.parse_energy(
        config.output_dir / f"{config.fragment1_name}_{index + 1:04d}.{config.output_extension}"
    )
    fragment2_scf = adapter.parse_energy(
        config.output_dir / f"{config.fragment2_name}_{index + 1:04d}.{config.output_extension}"
    )

    distances = [distance(frame[a - 1], frame[b - 1]) for a, b in config.distances]
    angles = [angle(frame[a - 1], frame[b - 1], frame[c - 1]) for a, b, c in config.angles]
    dihedrals = [dihedral(frame[a - 1], frame[b - 1], frame[c - 1], frame[d - 1]) for a, b, c, d in config.dihedrals]

    distortion_fragment1 = fragment1_scf - config.fragment1_energy
    distortion_fragment2 = fragment2_scf - config.fragment2_energy
    interaction_energy = complex_scf - (fragment1_scf + fragment2_scf)
    total_energy = complex_scf - config.fragment1_energy - config.fragment2_energy

    return AnalysisRow(
        step=index + 1,
        distances=distances,
        angles=angles,
        dihedrals=dihedrals,
        total_energy=total_energy * HARTREE_TO_KCAL,
        interaction_energy=interaction_energy * HARTREE_TO_KCAL,
        distortion_total=(distortion_fragment1 + distortion_fragment2) * HARTREE_TO_KCAL,
        distortion_fragment1=distortion_fragment1 * HARTREE_TO_KCAL,
        distortion_fragment2=distortion_fragment2 * HARTREE_TO_KCAL,
        complex_scf=complex_scf,
        fragment1_scf=fragment1_scf,
        fragment2_scf=fragment2_scf,
    )


def _append_analysis_row(config: ProjectConfig, row: AnalysisRow, initialize: bool) -> None:
    if initialize:
        headers = ["Step"]
        headers.extend(f"Distance {a}-{b}" for a, b in config.distances)
        headers.extend(f"Angle {a}-{b}-{c}" for a, b, c in config.angles)
        headers.extend(f"Dihedral {a}-{b}-{c}-{d}" for a, b, c, d in config.dihedrals)
        headers.extend(
            [
                "Total energy (kcal/mol)",
                "E(Int) (kcal/mol)",
                "E(Dist total) (kcal/mol)",
                f"E(Dist {config.fragment1_name}) (kcal/mol)",
                f"E(Dist {config.fragment2_name}) (kcal/mol)",
                "E(SCF complex) (hartree)",
                f"E(SCF {config.fragment1_name}) (hartree)",
                f"E(SCF {config.fragment2_name}) (hartree)",
            ]
        )
        config.analysis_path.write_text("\t".join(headers) + "\n", encoding="utf-8")

    values: list[str] = [str(row.step)]
    values.extend(f"{value:.5f}" for value in row.distances)
    values.extend(f"{value:.3f}" for value in row.angles)
    values.extend(f"{value:.3f}" for value in row.dihedrals)
    values.extend(
        [
            f"{row.total_energy:.5f}",
            f"{row.interaction_energy:.5f}",
            f"{row.distortion_total:.5f}",
            f"{row.distortion_fragment1:.5f}",
            f"{row.distortion_fragment2:.5f}",
            f"{row.complex_scf:.9f}",
            f"{row.fragment1_scf:.9f}",
            f"{row.fragment2_scf:.9f}",
        ]
    )
    with config.analysis_path.open("a", encoding="utf-8") as handle:
        handle.write("\t".join(values) + "\n")


def _cleanup(config: ProjectConfig) -> None:
    if not config.keep_input and config.input_dir.exists():
        shutil.rmtree(config.input_dir)
    if not config.keep_output and config.output_dir.exists():
        shutil.rmtree(config.output_dir)
    if not config.keep_log and config.log_path.exists():
        config.log_path.unlink()
