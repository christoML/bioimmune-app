import React, { useState } from "react";
import "../styles/bio.css";

const DataPreprocessing: React.FC = () => {
  const [open, setOpen] = useState(false);

  return (
    <div className="dp-container">
      {/* Main button */}
      <button
        className="bio-btn dp-main-btn"
        onClick={() => setOpen(!open)}
      >
        🧪 Data Preprocessing ▾
      </button>

      {open && (
        <div className="dp-dropdown">
          {/* Normalization */}
          <div className="dp-item">
            <span>.CEL Normalization</span>
            <span className="dp-arrow">▶</span>

            <div className="dp-submenu">
              <div className="dp-subitem">RMA</div>
            </div>
          </div>

          {/* Batch Effect */}
          <div className="dp-item">
            <span>Batch Effect Removal</span>
            <span className="dp-arrow">▶</span>

            <div className="dp-submenu">
              <div className="dp-subitem">ComBat</div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default DataPreprocessing;
