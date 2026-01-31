import React, { useState, useRef, useEffect, forwardRef, useImperativeHandle } from "react";
import "../styles/bio.css";

export interface GetMetadataRef {
  reset: () => void;
}

const GetMetadata = forwardRef<GetMetadataRef>((props, ref) => {
  const [logs, setLogs] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const logRef = useRef<HTMLDivElement>(null);

  useImperativeHandle(ref, () => ({
    reset: () => {
      setLogs([]);
      setLoading(false);
    },
  }));

  const handleGetMetadata = async () => {
    setLogs(["📦 Loading metadata from datasets..."]);
    setLoading(true);

    try {
      const res = await fetch("http://localhost:8000/get_geo_groups", {
        method: "POST",
      });

      if (!res.ok) throw new Error(`Backend returned ${res.status}`);
      const data = await res.json();

      if (data.success) {
        setLogs([...data.logs, "✅ Metadata processing completed!"]);
      } else {
        setLogs([`❌ Error: ${data.error}`]);
      }
    } catch (err: any) {
      setLogs([
        "❌ Failed to contact backend. Make sure the server is running.",
      ]);
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight;
  }, [logs]);

  return (
    <div className="gse-input-container">
      <button
        onClick={handleGetMetadata}
        className={`gse-input-button ${loading ? "downloading" : ""}`}
        disabled={loading}
      >
        {loading && <span className="spinner"></span>}
        {loading ? "Processing..." : "Get Metadata"}
      </button>

      <div ref={logRef} className="gse-input-logs">
        {logs.map((log, i) => (
          <div key={i}>{log}</div>
        ))}
      </div>
    </div>
  );
});

export default GetMetadata;
