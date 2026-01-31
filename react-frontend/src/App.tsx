import { BrowserRouter, Routes, Route, Navigate } from "react-router-dom";
import Toolbar from "./components/Toolbar";

import BioProcessing from "./pages/BioProcessing";

function App() {
  return (
    <BrowserRouter>
      <Toolbar />
      <Routes>
        <Route path="/" element={<Navigate to="/bio-processing" />} />
        <Route path="/bio-processing" element={<BioProcessing />} />
      </Routes>
    </BrowserRouter>
  );
}

export default App;
