from autodias.adapters.base import EngineAdapter
from autodias.adapters.gaussian import GaussianAdapter
from autodias.adapters.orca import OrcaAdapter


def get_adapter(name: str) -> EngineAdapter:
    normalized = name.lower()
    if normalized == "gaussian":
        return GaussianAdapter()
    if normalized == "orca":
        return OrcaAdapter()
    raise ValueError(f"Unsupported engine: {name}")
