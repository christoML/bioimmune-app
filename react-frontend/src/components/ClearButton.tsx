import React from "react";
import "../styles/edc.css";

interface ClearButtonProps {
  onReset: () => void;
}

const ClearButton: React.FC<ClearButtonProps> = ({ onReset }) => {
  return (
    <button className="edc-btn edc-btn-clear" onClick={onReset}>
      🔄 Clear All
    </button>
  );
};

export default ClearButton;
