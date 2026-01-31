import React, { useState, useRef, useEffect, forwardRef, useImperativeHandle } from "react";
import "../styles/bio.css";

export interface GSEInputRef {
  reset: () => void;
}

const GSEInput = forwardRef<GSEInputRef>((props, ref) => {
  const [input, setInput] = useState("");
  const [error, setError] = useState("");
  const [logs, setLogs] = useState<string[]>([]);
  const [downloading, setDownloading] = useState(false);
  const logRef = useRef<HTMLDivElement>(null);

  const gseRegex = /^(\s*GSE\d+\s*)(,\s*GSE\d+\s*)*$/i;

  useImperativeHandle(ref, () => ({
    reset: () => {
      setInput("");
      setError("");
      setLogs([]);
      setDownloading(false);
    },
  }));

  const handleSubmit = async () => {
    const trimmedInput = input.trim();
    if (!gseRegex.test(trimmedInput)) {
      setError(
        "Invalid format! Enter GSE IDs separated by commas, e.g., GSE12345,GSE67890"
      );
      setLogs([]);
      return;
    }

    setError("");
    setLogs(["📦 Download started, please wait..."]);
    setDownloading(true);

    const gseList = trimmedInput.split(",").map((id) => id.trim());

    try {
      const res = await fetch("http://localhost:8000/datasets", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ gse_ids: gseList.join(",") }),
      });

      if (!res.ok) throw new Error(`Backend returned ${res.status}`);

      const data = await res.json();

      if (data.success) {
        setLogs((prev) => [...prev, ...data.logs]);
        setLogs((prev) => [
          ...prev,
          "✅ GSE files download finished! Check the 'datasets' folder.",
        ]);
      } else {
        setError("Backend failed to download datasets.");
        setLogs([]);
      }
    } catch (err: any) {
      console.error(err);
      setError("Failed to contact backend. Make sure the server is running.");
      setLogs([]);
    } finally {
      setDownloading(false);
    }
  };

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight;
  }, [logs]);

  return (
    <div className="gse-input-container">
      <label htmlFor="gse-input" className="gse-input-label">
        Enter GSE IDs (comma-separated):
      </label>
      <input
        type="text"
        id="gse-input"
        value={input}
        onChange={(e) => setInput(e.target.value)}
        placeholder="GSE12345,GSE67890"
        className="gse-input-textbox"
        disabled={downloading}
      />
      <button
        onClick={handleSubmit}
        className={`gse-input-button ${downloading ? "downloading" : ""}`}
        disabled={downloading}
      >
        {downloading && <span className="spinner"></span>}
        {downloading ? "Downloading..." : "Submit"}
      </button>

      {error && <p className="gse-input-error">{error}</p>}

      <div ref={logRef} className="gse-input-logs">
        {logs.map((log, i) => (
          <div key={i}>{log}</div>
        ))}
      </div>
    </div>
  );
});

export default GSEInput;
