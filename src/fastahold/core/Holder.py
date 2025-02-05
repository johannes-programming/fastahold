from typing import *

from datahold import OkayList

from fastahold._utils import splitfiletext
from fastahold.core.Record import Record

__all__ = ["Holder"]


class Holder(OkayList):
    @property
    def data(self) -> Any:
        "This property represents the records."
        return list(self._data)

    @data.setter
    def data(self, value: Any) -> None:
        self._data = [Record(x) for x in value]

    @data.deleter
    def data(self) -> None:
        self._data = list()

    def dump(self, stream: "BinaryIO") -> None:
        "This method dumps the current instance into a wb stream."
        for rec in self:
            stream.write(rec.text.encode())

    def dumpintofile(self, file: Any) -> None:
        "This method dumps the current instance into a file."
        with open(file, "wb") as stream:
            self.dump(stream)

    def dumps(self) -> str:
        "This method dumps the current instance as a string."
        ans = ""
        for rec in self:
            ans += rec.text
        return ans

    @classmethod
    def load(cls, stream: "BinaryIO") -> Self:
        "This classmethod loads a new instance from an rb stream."
        text = stream.read().decode()
        ans = cls.loads(text)
        return ans

    @classmethod
    def loadfromfile(cls, file: Any) -> Self:
        "This classmethod loads a new instance from a file."
        with open(file, "rb") as stream:
            return cls.load(stream)

    @classmethod
    def loads(cls, string: Any) -> Self:
        "This classmethod loads a new instance from a string."
        code = splitfiletext(string)
        data = list()
        for block in code:
            description, seq = (block + "\n").split("\n", 1)
            rec = Record(description=description, seq=seq)
            data.append(rec)
        ans = cls(data)
        return ans
