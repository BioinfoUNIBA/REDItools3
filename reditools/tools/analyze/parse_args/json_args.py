"""Save and load CLI options using JSON files."""
import argparse
import json
from pathlib import Path

json_args_filename = "cli_args.json"

def args_to_json(
    options: argparse.Namespace,
    dirname: str,
    filename: str=json_args_filename,
) -> None:
    """Save commandline arguments to a JSON file.

    Parameters
    ----------
    options : argparse.Namespace
        The parsed commandline options.
    dirname : str
        Path to save file to.
    filename : str
        Name of the file (defaults to json_args_filename)
    """
    with Path(dirname, filename).open("w") as stream:
        json.dump(vars(options), stream)  # noqa: WPS421

def args_from_json(
    dirname: str,
    filename: str=json_args_filename,
) -> argparse.Namespace:
    """Load commandline arguments from a JSON file.

    Parameters
    ----------
    dirname : str
        Path to the JSON file folder
    filename : str
        JSON file to load arguments from (defaults to json_args_filename)

    Returns
    -------
    argparse.Namespace
        Commandline arguments for reditools analyze
    """
    with Path(dirname, filename).open("r") as stream:
        json_args = json.load(stream)
    return argparse.Namespace(**json_args)
