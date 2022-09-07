from dataclasses import dataclass
from typing import Optional
import pandas as pd

from .data_processor import DataSchema

@dataclass
class DataSource:
    _data: pd.DataFrame

    def create_pivot_table(self) -> pd.DataFrame:
        pt = self._data.pivot_table(
            values=DataSchema.AMOUNT,
            index=[DataSchema.CATEGORY],
            aggfunc="sum",
            fill_value=0,
            dropna=False,
        )
        return pt.reset_index().sort_values(DataSchema.AMOUNT, ascending=False)

    @property
    def row_count(self) -> int:
        return self._data.shape[0]

    @property
    def all_years(self) -> list[str]:
        return self._data[DataSchema.YEAR].tolist()

    @property
    def all_months(self) -> list[str]:
        return self._data[DataSchema.MONTH].tolist()

    @property
    def all_categories(self) -> list[str]:
        return self._data[DataSchema.CATEGORY].tolist()


