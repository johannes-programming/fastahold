from typing import *

from datahold import OkayList

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
