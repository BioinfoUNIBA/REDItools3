from __future__ import annotations

import argparse
import tempfile
from typing import Callable, Any

from reditools import reditools
from reditools.tools.analyze.parse_args.json_args import args_from_json

class ValueBelowMinimumError(argparse.ArgumentTypeError):
    def __init__(self, min_value: int | float) -> None:
        self.message = f"Value must be at least {min_value}."
        super().__init__(self.message)

class ValueAboveMaximumError(argparse.ArgumentTypeError):
    def __init__(self, max_value: int | float) -> None:
        self.message = f"Value cannot be larger than {max_value}."
        super().__init__(self.message)

class CastValueError(argparse.ArgumentTypeError):
    def __init__(self, typ: str, cli_val: Any) -> None:
        self.message = f"Invalid {typ} value: {cli_val}"
        super().__init__(self.message)

def check_number_bounds(
        number: float,
        min_value: float | None = None,
        max_value: float | None = None,
) -> None:
    """Check if a number is within specified bounds.

    Parameters
    ----------
    number : float
        The number to check.
    min_value : float | None, optional
        The minimum allowed value, by default None.
    max_value : float | None, optional
        The maximum allowed value, by default None.

    Raises
    ------
    ValueBelowMinimumError, ValueAboveMaximumError
        If the number is outside the specified bounds.
    """
    if min_value is not None and number < min_value:
        raise ValueBelowMinimumError(min_value)
    if max_value is not None and number > max_value:
        raise ValueAboveMaximumError(max_value)

def bounded_int(
        min_value: int | None = None,
        max_value: int | None = None,
) -> Callable:
    """Create a function that parses a string to a bounded integer.

    Parameters
    ----------
    min_value : int | None, optional
        The minimum allowed value, by default None.
    max_value : int | None, optional
        The maximum allowed value, by default None.

    Returns
    -------
    Callable
        A function that takes a string and returns a bounded integer.
    """
    def subfn(cli_value: str) -> int:  # noqa: WPS430
        try:
            int_value = int(cli_value)
        except ValueError as exc:
            raise CastValueError("int", cli_value) from exc
        check_number_bounds(int_value, min_value, max_value)
        return int_value
    return subfn


def bounded_float(
        min_value: float | None = None,
        max_value: float | None = None,
) -> Callable:
    """Create a function that parses a string to a bounded float.

    Parameters
    ----------
    min_value : float | None, optional
        The minimum allowed value, by default None.
    max_value : float | None, optional
        The maximum allowed value, by default None.

    Returns
    -------
    Callable
        A function that takes a string and returns a bounded float.
    """
    def subfn(cli_value: str) -> float:  # noqa: WPS430
        try:
            float_value = float(cli_value)
        except ValueError as exc:
            raise CastValueError("float", cli_value) from exc
        check_number_bounds(float_value, min_value, max_value)
        return float_value
    return subfn
