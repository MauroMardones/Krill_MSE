# Short-Cut MSE for Antarctic Krill – WAP

This project implements a **short-cut Management Strategy Evaluation (MSE)** to test the robustness of the **fixed TAC rule (145,000 t/year)** currently used by **CCAMLR** for Antarctic krill (*Euphausia superba*) in the Western Antarctic Peninsula (WAP).

The approach is adapted from the work of **Henning Winker (GFCM)** on **Blackspot seabream**, specifically his framework for shortcut MSEs in sex-structured models.

---

## Overview

- Simulates **true stock dynamics** using an operating model (OM).
- Emulates **assessment outcomes** ($SSB$, $F$) with error to derive TAC via the HCR.
- Avoids yearly re-fitting of the assessment model, reducing complexity and computation.

---

## Files

- `ShortCut_MSE_Krill.Rmd`: Main simulation script.
- `base.model1.4.rds`: Stock assessment output used to initialize the OM.
- `MSE_Krill_FLR.Rmd`: Additional modeling and testing code.
- `.gitignore`, `LICENSE`, `README.md`: Project setup and metadata.

---

## Tools

- R and the **[FLR framework](https://flr-project.org/)**

---

Let me know if you want a Spanish version too or add setup instructions!