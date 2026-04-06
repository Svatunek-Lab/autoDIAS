from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


IntTuple2 = tuple[int, int]
IntTuple3 = tuple[int, int, int]
IntTuple4 = tuple[int, int, int, int]


@dataclass(slots=True)
class ProjectConfig:
    name: str
    input_path: Path
    analysis_only: bool
    charge: int
    multiplicity: int
    reorder: bool
    reduce: bool
    reduce_threshold: float
    auto_detect_fragments: bool
    fragment1_name: str
    fragment2_name: str
    fragment1_charge: int
    fragment1_multiplicity: int
    fragment1_energy: float
    fragment2_charge: int
    fragment2_multiplicity: int
    fragment2_energy: float
    fragment1_atoms: list[int]
    fragment2_atoms: list[int]
    analysis_enabled: bool
    distances: list[IntTuple2]
    angles: list[IntTuple3]
    dihedrals: list[IntTuple4]
    keep_xyz: bool
    keep_input: bool
    keep_output: bool
    keep_log: bool
    prepare_only: bool
    engine_name: str
    input_extension: str
    output_extension: str
    command_template: str
    input_template: str
    base_dir: Path

    @property
    def input_dir(self) -> Path:
        return self.base_dir / f"{self.name}_input"

    @property
    def output_dir(self) -> Path:
        return self.base_dir / f"{self.name}_output"

    @property
    def xyz_dir(self) -> Path:
        return self.base_dir / f"{self.name}_xyz"

    @property
    def log_path(self) -> Path:
        suffix = "_analysisonly_log.txt" if self.analysis_only else "_log.txt"
        return self.base_dir / f"{self.name}{suffix}"

    @property
    def analysis_path(self) -> Path:
        return self.base_dir / f"{self.name}_DIAS.txt"


@dataclass(slots=True)
class StructureSet:
    atoms: list[str]
    frames: list[np.ndarray]
    titles: list[str]
    fragment1_atoms: list[int] = field(default_factory=list)
    fragment2_atoms: list[int] = field(default_factory=list)


@dataclass(slots=True)
class AnalysisRow:
    step: int
    distances: list[float]
    angles: list[float]
    dihedrals: list[float]
    total_energy: float
    interaction_energy: float
    distortion_total: float
    distortion_fragment1: float
    distortion_fragment2: float
    complex_scf: float
    fragment1_scf: float
    fragment2_scf: float
