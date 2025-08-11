from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal, TypeAlias

    import numpy as np

    from vasp_manager.calculation_manager.bulkmod import BulkmodCalculationManager
    from vasp_manager.calculation_manager.elastic import ElasticCalculationManager
    from vasp_manager.calculation_manager.rlx import RlxCalculationManager
    from vasp_manager.calculation_manager.rlx_coarse import RlxCoarseCalculationManager
    from vasp_manager.calculation_manager.static import StaticCalculationManager

    Filepath = str | Path
    Filepaths = list[Path]
    SourceDirectory = str | Path
    WorkingDirectory = str | Path
    Floating = float | np.floating
    CalculationType = Literal["rlx-coarse", "rlx", "static", "bulkmod", "elastic"]
    CalculationManager: TypeAlias = (
        BulkmodCalculationManager
        | ElasticCalculationManager
        | RlxCalculationManager
        | RlxCoarseCalculationManager
        | StaticCalculationManager
    )
