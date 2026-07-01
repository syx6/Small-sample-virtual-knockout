"""Prior-constrained virtual knockout toolkit."""

__all__ = ["VirtualKOResult", "run_virtual_ko"]


def __getattr__(name: str):
    if name in __all__:
        from .core import VirtualKOResult, run_virtual_ko

        return {"VirtualKOResult": VirtualKOResult, "run_virtual_ko": run_virtual_ko}[name]
    raise AttributeError(f"module 'vkx' has no attribute {name!r}")
