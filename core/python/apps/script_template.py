#!/usr/bin/env python
"""
Template and Example for a simple python app with argument parser.

https://docs.python.org/3.4/howto/argparse.html

"""


def main():
    """Entry function."""
    import argparse

    # intialisize parser
    parser = argparse.ArgumentParser(description="sample application")

    # define some options (better import from defaultParser)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Be verbose",
    )
    parser.add_argument(
        "-i",
        "--integer",
        dest="integer",
        type=int,
        help="Eine Integer Variable mit long und short option",
        default=0,
    )
    parser.add_argument(
        "-s",
        dest="string",
        type=str,
        help="Eine string Variable mit short option",
        default="none",
    )

    parser.add_argument("target")

    args = parser.parse_args()

    # results are in args
    if args.verbose:
        print(args)

    print("Integer is:", args.integer)
    print("String is:", args.string)

    print("Target is:", args.target)


# def main(...)

if __name__ == "__main__":
    main()
