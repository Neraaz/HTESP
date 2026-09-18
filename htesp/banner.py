#!/usr/bin/env python
"""The HTESP banner and citation, written to ``log`` by every command."""
from __future__ import annotations

__version__ = "2.0.0"

LOGO = r"""
                        ****    ****   ***********   *********   ********       ********
                        |  |    |  |       | |       | |____     | |            |  |  \ \
                        |  |____|  |       | |       |  ____|    | |*****       |  |__| |
                        |   ____   |       | |       | |_____          | |      |  _____/
                        |  |    |  |       | |       |_______|   ******| |      |_ |
                        |__|    |__|       |_| *********************************************
                        **********************
                                      High Throughput Electron-Structure Package
"""

AUTHORS = """
                                              Program written by

                                Niraj K Nepal, PhD       &       Lin-Lin Wang, PhD
                          Email: nnepal@ameslab.gov               llw@ameslab.gov
"""

#: the published reference.  The banner used to cite the arXiv preprint while
#: README.md cited the journal version; both now point at the same paper.
CITATION = (
    "N. K. Nepal, P. C. Canfield, and L.-L. Wang, HTESP (High-Throughput Electronic\n"
    "Structure Package): A package for high-throughput ab initio calculations,\n"
    "Comput. Mater. Sci. 244, 113247 (2024).  doi:10.1016/j.commatsci.2024.113247\n"
)

RULE = "_" * 131


def banner(version: str = __version__) -> str:
    """The full banner as one string."""
    return (
        f"{RULE}\n\n{LOGO}\n{AUTHORS}\n"
        f"                                            version {version}\n\n"
        "To support development, please cite the following paper and the papers\n"
        "referenced therein for the calculations you run.\n\n"
        f"{CITATION}\n{RULE}\n"
    )
