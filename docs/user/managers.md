# Managers

VaspManager handles the creation of a set of CalculationManagers for each material
and specified calculation types.
```mermaid
graph TD
A[VaspManager] --> C(Material 1) & D(Material 2) & E(Material ...)
subgraph identifier[" "]
direction LR
  F(Material N) --> G[RlxCoarseManager]
  F --> H[RlxManager]
  F --> I[StaticManager]
  F --> J[BulkmodManager]
  F --> K[ElasticManager]
end
A[VaspManager] --> identifier[" "]
```

Each CalculationManager has its own logical flow in order to handle VASP input creation,
job submission and monitoring, and result processing.
```mermaid
graph TD
A[CalculationManager] --> B[VaspInputCreator]
A --> C[JobManager]
A -.-> D[Analzyer]
```
