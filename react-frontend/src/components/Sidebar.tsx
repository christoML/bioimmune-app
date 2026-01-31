import React, { useRef } from "react";
import "../styles/edc.css";
import pliadesLogo from "../assets/pliades.png";

// Action components (unchanged behavior)
import CreateAsset from "../components/CreateAsset";
import CreatePolicy from "../components/CreatePolicy";
import CreateContractDefinition from "../components/CreateContractDefinition";
import FetchCatalog from "../components/FetchCatalog";
import NegotiateContract from "../components/NegotiateContract";
import GetContractNegotiation from "../components/GetContractNegotiation";
import TransferData from "../components/TransferData";

const Sidebar: React.FC = () => {
  const createAssetRef = useRef<any>(null);
  const createPolicyRef = useRef<any>(null);
  const createContractRef = useRef<any>(null);
  const fetchCatalogRef = useRef<any>(null);
  const negotiateContractRef = useRef<any>(null);
  const getContractNegotiationRef = useRef<any>(null);
  const transferDataRef = useRef<any>(null);

  return (
    <div className="sidebar">
      {/* Logo (unchanged) */}
      <div className="sidebar-logo">
        <img src={pliadesLogo} alt="PLIADES" />
      </div>

      <div className="sidebar-menu">
        {/* Consume */}
        <div className="sidebar-section">Consume</div>
        <FetchCatalog ref={fetchCatalogRef} />

        {/* Provide */}
        <div className="sidebar-section">Provide</div>
        <CreateAsset ref={createAssetRef} />
        <CreatePolicy ref={createPolicyRef} />
        <CreateContractDefinition ref={createContractRef} />
        <NegotiateContract ref={negotiateContractRef} />
        <GetContractNegotiation ref={getContractNegotiationRef} />
        <TransferData ref={transferDataRef} />
      </div>
    </div>
  );
};

export default Sidebar;
