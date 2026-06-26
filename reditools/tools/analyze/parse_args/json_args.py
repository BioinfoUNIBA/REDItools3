import argparse
import json
import os

json_args_filename = 'cli_args.json'

def args_to_json(
    options: argparse.Namespace,
    dirname: str,
    filename: str=json_args_filename,
) -> None:
    """
    Save commandline arguments to a JSON file.

    Parameters:
        options : argparse.Namespace
            The parsed commandline options.
        filename : str
            Path to save arguments to.
    """
    with open(os.path.join(dirname, filename), 'w') as stream:
        json.dump(vars(options), stream)  # noqa: WPS421

def args_from_json(
    dirname: str,
    filename: str=json_args_filename,
) -> argparse.Namespace:
    """
    Load commandline arguments from a JSON file.

    Parameters
    ----------
    filename : str
        JSON file to load arguments from

    Returns
    -------
    argparse.Namespace
        Commandline arguments for reditools analyze
    """
    with open(os.path.join(dirname, filename), 'r') as stream:
        json_args = json.load(stream)
    return argparse.Namespace(**json_args)
