import biologo from "../assets/bioinformatics_logo.png";

const LogoBanner = () => {
  return (
    <div style={styles.container}>
      <img
        src={biologo}
        alt="BioImmune Logo"
        style={styles.image}
      />
    </div>
  );
};

const styles = {
  container: {
    marginTop: "60px",
    textAlign: "center" as const,
  },
  image: {
    width: "240px",
    marginBottom: "16px",
  },
  subtitle: {
    color: "#3a86ff",
    fontWeight: 300,
  },
};

export default LogoBanner;
