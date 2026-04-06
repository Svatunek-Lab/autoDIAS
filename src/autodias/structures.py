from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from autodias.geometry import rmsd_after_alignment
from autodias.models import ProjectConfig, StructureSet


GAUSSIAN_COORD_RE = re.compile(
    r"^\s*\d+\s+(\d+)\s+\d+\s+([-0-9.Ee]+)\s+([-0-9.Ee]+)\s+([-0-9.Ee]+)\s*$"
)


def load_structures(config: ProjectConfig) -> StructureSet:
    suffix = config.input_path.suffix.lower()
    if suffix == ".xyz":
        structures = _read_xyz(config.input_path)
    else:
        structures = _read_gaussian_output(config.input_path)

    if config.fragment1_atoms and config.fragment2_atoms:
        structures.fragment1_atoms = sorted(config.fragment1_atoms)
        structures.fragment2_atoms = sorted(config.fragment2_atoms)
        _validate_fragments(structures)
    elif config.auto_detect_fragments:
        frag1, frag2 = detect_fragments(structures, config.charge)
        structures.fragment1_atoms = frag1
        structures.fragment2_atoms = frag2
    else:
        raise ValueError("Fragments must be supplied if auto detection is disabled")

    if config.reorder:
        reorder_frames(structures)
    if config.reduce:
        reduce_frames(structures, config.reduce_threshold)

    return structures


def _read_xyz(path: Path) -> StructureSet:
    lines = path.read_text().splitlines()
    idx = 0
    atoms: list[str] = []
    frames: list[np.ndarray] = []
    titles: list[str] = []
    while idx < len(lines):
        if not lines[idx].strip():
            idx += 1
            continue
        natoms = int(lines[idx].strip())
        title = lines[idx + 1].strip()
        block = lines[idx + 2 : idx + 2 + natoms]
        frame_atoms: list[str] = []
        coords: list[list[float]] = []
        for row in block:
            parts = row.split()
            frame_atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if not atoms:
            atoms = frame_atoms
        elif atoms != frame_atoms:
            raise ValueError("Atom ordering changes between XYZ frames")
        frames.append(np.asarray(coords, dtype=float))
        titles.append(title or str(len(frames)))
        idx += natoms + 2
    if not frames:
        raise ValueError(f"No XYZ frames found in {path}")
    return StructureSet(atoms=atoms, frames=frames, titles=titles)


def _read_gaussian_output(path: Path) -> StructureSet:
    atom_numbers: list[int] | None = None
    atoms: list[str] = []
    frames: list[np.ndarray] = []
    titles: list[str] = []
    current_block: list[str] = []
    in_orientation = False
    dash_count = 0

    for line in path.read_text(errors="replace").splitlines():
        if "Standard orientation:" in line or "Input orientation:" in line:
            in_orientation = True
            current_block = []
            dash_count = 0
            continue
        if not in_orientation:
            continue
        if "-----" in line:
            dash_count += 1
            if dash_count == 2:
                current_block = []
            elif dash_count == 3:
                frame, atom_numbers = _parse_gaussian_block(current_block, atom_numbers)
                if atom_numbers and not atoms:
                    atoms = [Chem.GetPeriodicTable().GetElementSymbol(n) for n in atom_numbers]
                frames.append(frame)
                titles.append(str(len(frames)))
                in_orientation = False
            continue
        if dash_count >= 2:
            current_block.append(line)

    if not frames:
        raise ValueError(f"No Gaussian coordinate blocks found in {path}")
    return StructureSet(atoms=atoms, frames=frames, titles=titles)


def _parse_gaussian_block(lines: list[str], atom_numbers: list[int] | None) -> tuple[np.ndarray, list[int]]:
    coords: list[list[float]] = []
    parsed_atom_numbers: list[int] = [] if atom_numbers is None else atom_numbers[:]
    for line in lines:
        match = GAUSSIAN_COORD_RE.match(line)
        if not match:
            continue
        atomic_number = int(match.group(1))
        x = float(match.group(2))
        y = float(match.group(3))
        z = float(match.group(4))
        coords.append([x, y, z])
        if atom_numbers is None:
            parsed_atom_numbers.append(atomic_number)
    if not coords:
        raise ValueError("Encountered empty Gaussian coordinate block")
    return np.asarray(coords, dtype=float), parsed_atom_numbers


def detect_fragments(structures: StructureSet, charge: int) -> tuple[list[int], list[int]]:
    fragment_splits: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for frame in structures.frames:
        components = _connected_components(structures.atoms, frame, charge)
        if len(components) < 2:
            continue
        if len(components) == 2:
            ordered = tuple(sorted(tuple(sorted(comp)) for comp in components))
            fragment_splits.append((ordered[0], ordered[1]))
            continue
        first = tuple(sorted(components[0]))
        second = tuple(sorted(idx for comp in components[1:] for idx in comp))
        fragment_splits.append((first, second))

    if not fragment_splits:
        raise ValueError("Could not determine fragments automatically from the available structures")

    common = Counter(fragment_splits).most_common(1)[0][0]
    return list(common[0]), list(common[1])


def _connected_components(atoms: list[str], frame: np.ndarray, charge: int) -> list[list[int]]:
    xyz_lines = [str(len(atoms)), "generated by autoDIAS 2.0"]
    xyz_lines.extend(
        f"{atom} {coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}"
        for atom, coords in zip(atoms, frame, strict=True)
    )
    mol = Chem.MolFromXYZBlock("\n".join(xyz_lines))
    if mol is None:
        raise ValueError("RDKit could not build a molecule from the structure")
    rw_mol = Chem.RWMol(mol)
    rdDetermineBonds.DetermineConnectivity(rw_mol, charge=charge)
    adjacency = Chem.GetAdjacencyMatrix(rw_mol)

    seen: set[int] = set()
    components: list[list[int]] = []
    for start in range(len(atoms)):
        if start in seen:
            continue
        stack = [start]
        component: list[int] = []
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            component.append(node + 1)
            neighbors = np.where(adjacency[node] > 0)[0]
            stack.extend(int(n) for n in neighbors if int(n) not in seen)
        components.append(sorted(component))
    return sorted(components, key=lambda item: (len(item), item))


def _validate_fragments(structures: StructureSet) -> None:
    natoms = len(structures.atoms)
    all_atoms = sorted(structures.fragment1_atoms + structures.fragment2_atoms)
    if all_atoms != list(range(1, natoms + 1)):
        raise ValueError("Fragment atom assignments must cover each atom exactly once")


def reorder_frames(structures: StructureSet) -> None:
    if len(structures.frames) <= 2:
        return
    rmsds = [
        rmsd_after_alignment(structures.frames[idx], structures.frames[idx + 1])
        for idx in range(len(structures.frames) - 1)
    ]
    mean_rmsd = float(np.mean(rmsds))
    if mean_rmsd == 0.0:
        return
    spike_idx = int(np.argmax(rmsds))
    if rmsds[spike_idx] <= 3.0 * mean_rmsd:
        return
    left = structures.frames[: spike_idx + 1]
    right = structures.frames[spike_idx + 1 :]
    if not right:
        return
    if rmsd_after_alignment(left[0], right[0]) < rmsd_after_alignment(left[-1], right[0]):
        structures.frames = list(reversed(left)) + right
        structures.titles = list(reversed(structures.titles[: spike_idx + 1])) + structures.titles[spike_idx + 1 :]


def reduce_frames(structures: StructureSet, threshold: float) -> None:
    if len(structures.frames) <= 1:
        return
    kept_frames = [structures.frames[0]]
    kept_titles = [structures.titles[0]]
    for title, frame in zip(structures.titles[1:], structures.frames[1:], strict=True):
        if rmsd_after_alignment(kept_frames[-1], frame) >= threshold:
            kept_frames.append(frame)
            kept_titles.append(title)
    structures.frames = kept_frames
    structures.titles = kept_titles
