import * as React from "react"; // Use the namespace import
import "../styles/bio.css"; // Keep this for styling
import LogoBanner from "../components/LogoBanner";
import GSEInput from "../components/GSEInput";
import DataPreprocessing from "../components/DataPreprocessing";


//remove the warning about React
const { useState } = React;

// Clear All button component
const ClearButton = ({ onReset }: { onReset: () => void }) => (
  <button className="bio-btn bio-btn-clear" onClick={onReset}>
    🔄 Clear All
  </button>
);

function BioProcessing() {
  const [terminatedAlert, setTerminatedAlert] = useState("");
  const [gseResetKey, setGseResetKey] = useState(0); // force GSEInput reset

  // Clear All action
  const handleClearAll = () => {
    // increment key to remount GSEInput
    setGseResetKey(prev => prev + 1);

    // Show temporary alert
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
          margin: "20px 0"
        }}
      >
        <ClearButton onReset={handleClearAll} />
        {terminatedAlert && (
          <div
            className="bio-alert bio-alert-success"
            style={{ marginLeft: 10 }}
          >
            {terminatedAlert}
          </div>
        )}
      </div>

      {/* GSE Input Component with reset key */}
      <GSEInput key={gseResetKey} />
      <DataPreprocessing />
    </div>
  );
}

export default BioProcessing;
