from typing import *

__all__ = ["parse"]


def splitfiletext(self, value: Any) -> List[str]:
    "This function splits the file text."
    code = str(value)
    code = "\n" + code + "\n"
    code = code.split("\n>")
    return code
