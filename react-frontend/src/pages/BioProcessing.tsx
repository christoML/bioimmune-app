import * as React from "react";
import "../styles/bio.css";
import LogoBanner from "../components/LogoBanner";
import GSEInput from "../components/GSEInput";
import GetMetadata from "../components/GetMetadata";
import type { GetMetadataRef } from "../components/GetMetadata";
import type { GSEInputRef } from "../components/GSEInput";
import DataPreprocessing from "../components/DataPreprocessing";

const { useState, useRef } = React;

const ClearButton = ({ onReset }: { onReset: () => void }) => (
  <button className="bio-btn bio-btn-clear" onClick={onReset}>
    🔄 Clear All
  </button>
);

function BioProcessing() {
  const [terminatedAlert, setTerminatedAlert] = useState("");

  // Refs for child components
  const gseInputRef = useRef<GSEInputRef>(null);
  const metadataRef = useRef<GetMetadataRef>(null);

  const handleClearAll = () => {
    gseInputRef.current?.reset();
    metadataRef.current?.reset();

    setTerminatedAlert("All inputs cleared!");
    setTimeout(() => setTerminatedAlert(""), 3000);
  };

  return (
    <div style={{ padding: "50px", textAlign: "center" }}>
      <h1>BioImmune App</h1>
      <LogoBanner />

      {/* Clear All button + alert */}
      <div
        style={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          margin: "20px 0",
        }}
      >
        <ClearButton onReset={handleClearAll} />
        {terminatedAlert && (
          <div className="bio-alert bio-alert-success" style={{ marginLeft: 10 }}>
            {terminatedAlert}
          </div>
        )}
      </div>

      {/* Step 1: GSE Input */}
      <GSEInput ref={gseInputRef} />

      {/* Step 2: Get Metadata */}
      <GetMetadata ref={metadataRef} />

      {/* Step 3: Data Preprocessing */}
      <DataPreprocessing />
    </div>
  );
}

export default BioProcessing;
