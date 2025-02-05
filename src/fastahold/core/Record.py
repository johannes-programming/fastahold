from typing import *

from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from datarepr import datarepr
from overloadable import overloadable

from fastahold._utils import splitfiletext

__all__ = ["Record"]


class Record:
    __slots__ = ("_description", "_seq")

    @overloadable
    def __init__(self, *args: Any, **kwargs: Any) -> str:
        if len(args) != 0:
            return "other"
        if len(kwargs) != 1:
            return "fields"
        (key,) = kwargs.keys()
        if key in ("text", "bio"):
            return key
        return "fields"

    @__init__.overload("other")
    def __init__(self, other: Self, /) -> None:
        if not isinstance(other, type(self)):
            raise TypeError
        self.description = other.description
        self.seq = other.seq

    @__init__.overload("fields")
    def __init__(self, *, description: Any = "", seq: Any = "") -> None:
        self.description = description
        self.seq = seq

    @__init__.overload("bio")
    def __init__(self, *, bio: SeqRecord) -> None:
        self.bio = bio

    @__init__.overload("text")
    def __init__(self, *, text: str) -> None:
        self.text = text

    def __repr__(self) -> str:
        "This magic method implements repr(self)."
        return datarepr(
            type(self).__name__,
            description=self.description,
            seq=self.seq,
        )

    @property
    def bio(self) -> SeqRecord:
        "This property holds a corresponding Bio.SeqRecord.SeqRecord."
        return SeqRecord(
            seq=self.seq,
            id=self.id,
            name=self.id,
            description=self.description,
        )

    @bio.setter
    def bio(self, value: SeqRecord) -> None:
        self.text = FastaIO.as_fasta(value)

    @bio.deleter
    def bio(self) -> None:
        del self.text

    def copy(self) -> Self:
        "This method returns a copy of the current instance."
        return type(self)(self)

    @property
    def description(self) -> str:
        "This property holds the description of the fasta record."
        return self._description

    @description.setter
    def description(self, value: Any) -> None:
        self._description = value

    @description.deleter
    def description(self) -> None:
        del self._description

    @property
    def id(self) -> str:
        "This property holds the id of the record."
        if self.description == "":
            return ""
        return self._description.split(None, 1)[0]

    @id.setter
    def id(self, value: Any) -> None:
        value = str(value)
        self.description = self.description.replace(self.id, value, 1)

    @id.deleter
    def id(self) -> None:
        self.id = ""

    @property
    def seq(self) -> Seq:
        "This property holds the seq of the record."
        return self._seq

    @seq.setter
    def seq(self, value: Any) -> None:
        seq = Seq(value)
        seq = seq.replace("\n", "").replace("\r", "").replace(" ", "")
        self._seq = seq

    @seq.deleter
    def seq(self) -> None:
        self._seq = Seq("")

    @property
    def supplementary(self) -> str:
        "This property holds the supplementary information of the record."
        if self.description == "":
            return ""
        return self.description.split(None, 1)[1]

    @supplementary.setter
    def supplementary(self, value: Any) -> None:
        value = str(value)
        index = len(self.description) - len(self.supplementary)
        self.description = self.description[:index] + value

    @supplementary.deleter
    def supplementary(self) -> None:
        self.supplementary = ""

    @property
    def text(self) -> str:
        "This property holds a text representation of the record."
        ans = ">%s\n" % self.description
        for i in range(0, len(self.seq), 60):
            ans += "%s\n" % self.seq[i : i + 60]
        return ans

    @text.setter
    def text(self, value: Any) -> None:
        code: list[str] = splitfiletext(value)
        code = code[-1].split("\n")
        code = [x.rstrip() for x in code]
        description0 = self.description
        self.description = code.pop(0)
        try:
            self.seq = "".join(code)
        except:
            self.description = description0
            raise

    @text.deleter
    def text(self) -> None:
        del self.description
        del self.seq
