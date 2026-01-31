import { useState } from "react";
import { useNavigate } from "react-router-dom";

interface ToolbarProps {
  title?: string;
}

const Toolbar = ({ title = "BioImmune Dashboard" }: ToolbarProps) => {
  const navigate = useNavigate();
  const [open, setOpen] = useState(false);

  return (
    <header style={styles.header}>
      <div style={styles.left}>
        <strong>{title}</strong>
      </div>

      <nav style={styles.right}>
        <button style={styles.button} onClick={() => navigate("/")}>
          Dashboard
        </button>

        {/* EDC Samples dropdown */}
        <div style={{ position: "relative" }}>
          <button
            style={styles.button}
            onClick={() => setOpen((prev) => !prev)}
          >
            BioImmune Menu ▾
          </button>

          {open && (
            <div style={styles.dropdown}>
              <div
                style={styles.dropdownItem}
                onClick={() => {
                  navigate("/edc-samples-file-transfer-cloud");
                  setOpen(false);
                }}
              >
                File Transfer Cloud
              </div>

              <div
                style={styles.dropdownItem}
                onClick={() => {
                  navigate("/scalability-test");
                  setOpen(false);
                }}
              >
                Scalability Test
              </div>

              <div
              style={styles.dropdownItem}
              onClick={() => {
                navigate("/connector/test");
                setOpen(false);
              }}
              >
                Connector UI
              </div>


            </div>
          )}
        </div>

        <button style={styles.button} onClick={() => navigate("/about")}>
          About
        </button>
      </nav>
    </header>
  );
};

const styles = {
  header: {
    height: "60px",
    display: "flex",
    justifyContent: "space-between",
    alignItems: "center",
    padding: "0 24px",
    backgroundColor: "#0b132b",
    color: "white",
  },
  left: {
    fontSize: "1.2rem",
  },
  right: {
    display: "flex",
    gap: "12px",
  },
  button: {
    background: "transparent",
    color: "white",
    border: "1px solid #5bc0be",
    padding: "6px 12px",
    cursor: "pointer",
  },

  /* 🔽 Dropdown styles (isolated, no changes to existing layout) */
  dropdown: {
    position: "absolute" as const,
    top: "38px",
    right: 0,
    backgroundColor: "#0b132b",
    border: "1px solid #5bc0be",
    minWidth: "200px",
    zIndex: 1000,
  },
  dropdownItem: {
    padding: "8px 12px",
    cursor: "pointer",
    borderBottom: "1px solid #5bc0be",
  },
};

export default Toolbar;
