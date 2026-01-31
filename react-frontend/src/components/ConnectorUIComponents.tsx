const ConnectorOverview = ({ connectorId }: { connectorId: string }) => (
  <>
    <h3>Dashboard</h3>
    <p>Overview for connector <strong>{connectorId}</strong></p>
    <ul>
      <li>Assets: (coming soon)</li>
      <li>Policies: (coming soon)</li>
      <li>Contracts: (coming soon)</li>
    </ul>
  </>
);

const ConnectorAssets = ({ connectorId }: { connectorId: string }) => (
  <>
    <h3>Assets</h3>
    <p>Manage assets for <strong>{connectorId}</strong></p>
    {/* Later: reuse CreateAsset component here */}
  </>
);

const ConnectorPolicies = ({ connectorId }: { connectorId: string }) => (
  <>
    <h3>Policies</h3>
    <p>Manage policies for <strong>{connectorId}</strong></p>
    {/* Later: reuse CreatePolicy */}
  </>
);

const ConnectorContracts = ({ connectorId }: { connectorId: string }) => (
  <>
    <h3>Contracts</h3>
    <p>Manage contracts for <strong>{connectorId}</strong></p>
    {/* Later: CreateContractDefinition */}
  </>
);

const ConnectorTransfers = ({ connectorId }: { connectorId: string }) => (
  <>
    <h3>Transfers</h3>
    <p>Data transfers for <strong>{connectorId}</strong></p>
    {/* Later: TransferData */}
  </>
);

export {
  ConnectorOverview,
  ConnectorAssets,
  ConnectorPolicies,
  ConnectorContracts,
  ConnectorTransfers,
};
