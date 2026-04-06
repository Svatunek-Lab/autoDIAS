from __future__ import annotations

import tomllib
from pathlib import Path

from autodias.models import ProjectConfig


def _tuple_list(values: list[list[int]], size: int) -> list[tuple[int, ...]]:
    tuples: list[tuple[int, ...]] = []
    for value in values:
        if len(value) != size:
            raise ValueError(f"Expected {size} integers, got {value!r}")
        tuples.append(tuple(int(item) for item in value))
    return tuples


def load_config(path: str | Path) -> ProjectConfig:
    config_path = Path(path).resolve()
    with config_path.open("rb") as handle:
        raw = tomllib.load(handle)

    project = raw["project"]
    structure = raw["structure"]
    fragments = raw["fragments"]
    analysis = raw["analysis"]
    files = raw["files"]
    engine = raw["engine"]

    return ProjectConfig(
        name=project["name"],
        input_path=(config_path.parent / project["input_path"]).resolve(),
        analysis_only=bool(project.get("analysis_only", False)),
        charge=int(structure["charge"]),
        multiplicity=int(structure["multiplicity"]),
        reorder=bool(structure.get("reorder", False)),
        reduce=bool(structure.get("reduce", False)),
        reduce_threshold=float(structure.get("reduce_threshold", 0.2)),
        auto_detect_fragments=bool(fragments.get("auto_detect", True)),
        fragment1_name=fragments.get("fragment1_name", "fragment1"),
        fragment2_name=fragments.get("fragment2_name", "fragment2"),
        fragment1_charge=int(fragments.get("fragment1_charge", 0)),
        fragment1_multiplicity=int(fragments.get("fragment1_multiplicity", 1)),
        fragment1_energy=float(fragments["fragment1_energy"]),
        fragment2_charge=int(fragments.get("fragment2_charge", 0)),
        fragment2_multiplicity=int(fragments.get("fragment2_multiplicity", 1)),
        fragment2_energy=float(fragments["fragment2_energy"]),
        fragment1_atoms=[int(x) for x in fragments.get("fragment1_atoms", [])],
        fragment2_atoms=[int(x) for x in fragments.get("fragment2_atoms", [])],
        analysis_enabled=bool(analysis.get("enabled", True)),
        distances=_tuple_list(analysis.get("distances", []), 2),
        angles=_tuple_list(analysis.get("angles", []), 3),
        dihedrals=_tuple_list(analysis.get("dihedrals", []), 4),
        keep_xyz=bool(files.get("keep_xyz", True)),
        keep_input=bool(files.get("keep_input", True)),
        keep_output=bool(files.get("keep_output", True)),
        keep_log=bool(files.get("keep_log", True)),
        prepare_only=bool(files.get("prepare_only", False)),
        engine_name=engine["name"].lower(),
        input_extension=engine.get("input_extension", "inp"),
        output_extension=engine.get("output_extension", "out"),
        command_template=engine["command"],
        input_template=engine["input_template"],
        base_dir=config_path.parent.resolve(),
    )
