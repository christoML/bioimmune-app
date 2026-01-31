import React, { useState, useRef, useEffect } from "react";
import "../styles/bio.css";

const DataPreprocessing: React.FC = () => {
  const [open, setOpen] = useState(false);
  const [logs, setLogs] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const logRef = useRef<HTMLDivElement>(null);

  const runRMA = async (type: "CEL" | "series_matrix") => {
    setLogs([`📦 Running RMA on ${type}...`]);
    setLoading(true);

    try {
      const res = await fetch(
        `http://localhost:8000/gene_mat_preprocess?type=${type}`,
        { method: "POST" }
      );

      if (!res.ok) throw new Error(`Backend returned ${res.status}`);
      const data = await res.json();

      if (data.success) {
        setLogs([...data.logs, "✅ Gene matrices generated!"]);
      } else {
        setLogs([`❌ Error: ${data.error}`]);
      }
    } catch (err) {
      console.error(err);
      setLogs(["❌ Failed to contact backend."]);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (logRef.current) {
      logRef.current.scrollTop = logRef.current.scrollHeight;
    }
  }, [logs]);

  return (
    <div className="dp-container">
      <button
        className="bio-btn dp-main-btn"
        onClick={() => setOpen(!open)}
      >
        🧪 Data Preprocessing ▾
      </button>

      {open && (
        <div className="dp-dropdown">

          {/* CEL Normalization */}
          <div className="dp-item">
            <span>.CEL Normalization</span>
            <span className="dp-arrow">▶</span>

            <div className="dp-submenu">

              {/* RMA level */}
              <div className="dp-item">
                <span>RMA</span>
                <span className="dp-arrow">▶</span>

                <div className="dp-submenu">

                  <div
                    className="dp-subitem"
                    onClick={() => !loading && runRMA("CEL")}
                    style={{ fontWeight: "bold", cursor: "pointer" }}
                  >
                    {loading ? "Processing..." : "RMA on .CEL files"}
                  </div>

                  <div
                    className="dp-subitem"
                    onClick={() => !loading && runRMA("series_matrix")}
                    style={{ fontWeight: "bold", cursor: "pointer" }}
                  >
                    {loading ? "Processing..." : "RMA on series_matrix"}
                  </div>

                </div>
              </div>

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

      {logs.length > 0 && (
        <div ref={logRef} className="gse-input-logs">
          {logs.map((log, i) => (
            <div key={i}>{log}</div>
          ))}
        </div>
      )}
    </div>
  );
};

export default DataPreprocessing;
