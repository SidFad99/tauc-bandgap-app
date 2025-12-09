UV–Vis–NIR Tauc Band Gap App

Streamlit app to convert UV–Vis–NIR data into absorption-related quantities and build Tauc plots for direct and indirect optical band gap estimation.

The app is designed for researchers and students working with glass, semiconductors, and other optical materials.

--------------------------------------------------
1. Features
--------------------------------------------------

- Accepts Absorbance (A), Transmittance (T), or Reflectance (R) data
- Converts to:
  - Absorption coefficient α (for A or T)
  - Kubelka–Munk remission function F(R) as a pseudo-absorption term (for diffuse reflectance)
- Computes:
  - Photon energy hν (eV) from wavelength (nm)
  - Tauc variables:
    - (αhν)² for allowed direct transitions
    - (αhν)¹ᐟ² for allowed indirect transitions
- Interactive Tauc plots with:
  - Sliders to choose the quasi-linear region
  - Automatic linear regression
  - Estimated Eg(direct) and Eg(indirect) from x-intercepts
  - R² values and short comments on fit quality
- Extra teaching tools:
  - “Method and equations” expander with the key relations in LaTeX
  - Automatically generated methods paragraph for lab reports or papers
- Export of processed data (α, αhν, Tauc columns) as CSV

--------------------------------------------------
2. Requirements
--------------------------------------------------

- Python 3.8 or newer

Python packages:

    pip install streamlit pandas numpy matplotlib

--------------------------------------------------
3. Installation
--------------------------------------------------

1. Clone the repository

    git clone https://github.com/<your-username>/tauc-bandgap-app.git
    cd tauc-bandgap-app

2. (Optional) create and activate a virtual environment

Using venv:

    python -m venv venv
    # Windows
    venv\Scripts\activate
    # macOS / Linux
    source venv/bin/activate

3. Install dependencies

If you have a requirements.txt:

    pip install -r requirements.txt

If not, run:

    pip install streamlit pandas numpy matplotlib

--------------------------------------------------
4. Running the app
--------------------------------------------------

From the repository folder:

    streamlit run tauc_bandgap_app_v4.py

Streamlit will open the app in your browser (usually at http://localhost:8501).

--------------------------------------------------
5. Input data format
--------------------------------------------------

The app expects tabular data in CSV or TXT form with at least:

- One column with wavelength in nm
- One column with Absorbance / Transmittance / Reflectance

Example files:

Wavelength,Abs
250,0.0123
251,0.0131
...

Wavelength,T_percent
250,72.3
251,71.8
...

Wavelength,R_percent
250,9.52
251,9.48
...

The first row should contain column names. In the app you then choose:

- Wavelength column (nm)
- Data column (Abs, T, or R)

--------------------------------------------------
6. Workflow inside the app
--------------------------------------------------

Step 0 – Measurement settings (sidebar)
- Choose data type: Absorbance A, Transmittance T, or Reflectance R (Kubelka–Munk)
- If A or T: enter sample thickness and units
- If T or R: choose whether the data are fraction (0–1) or percent (0–100)

Step 1 – Upload data
- Upload a CSV or TXT file exported from your spectrophotometer or spreadsheet software.

Step 2–4 – Check preview and calculated quantities
- Confirm that the preview table looks sensible.
- The app converts to photon energy hν, absorption-like quantity, αhν, (αhν)², and (αhν)¹ᐟ².
- A short data summary reports wavelength range, energy range, NaNs and negative values.

Step 5 – Inspect spectra and Tauc plots
- Use three tabs:
  - Raw spectrum (signal vs wavelength)
  - Tauc (direct): (αhν)² vs hν
  - Tauc (indirect): (αhν)¹ᐟ² vs hν

Step 6 – Choose energy windows and extract Eg
- For both direct and indirect fits:
  - Move the slider to select the straight part near the absorption edge.
  - The app performs a linear fit and shows:
    - Eg (x-intercept)
    - Fit equation
    - R² value
    - Short comment on linearity
  - Detailed plots show:
    - Blue curve = Tauc data
    - Orange line = linear fit
    - Dashed line = Eg

Step 7 – Summary and report text
- The app displays Direct Eg and Indirect Eg as Streamlit metrics.
- A ready-made methods paragraph is generated, including:
  - Measurement mode (A / T / R)
  - Thickness statement
  - Tauc exponents
  - Fit ranges
  - Eg values and R²

Step 8 – Download data
- Download a CSV with hv_eV, alpha_like, alpha_hv, tauc_direct, and tauc_indirect for further analysis.

--------------------------------------------------
7. Theory (short overview)
--------------------------------------------------

The app follows a standard approach for optical band gap determination:

- Photon energy and wavelength are related by

      hν [eV] = 1240 / λ [nm]

- For transmission data, the Beer–Lambert law connects absorbance A and absorption coefficient α:

      α = (2.303 A) / d

  where d is the sample thickness in cm.

- For diffuse reflectance, the Kubelka–Munk remission function

      F(R∞) = (1 − R∞)² / (2 R∞)

  is used as a pseudo-absorption term for strongly scattering samples.

- The Tauc relation describes the near-edge dependence of the absorption:

      (αhν)^n ∝ (hν − Eg)

  with exponent n = 2 for allowed direct transitions and n = 1/2 for allowed indirect transitions.

By plotting (αhν)^n versus hν and extrapolating the linear region to zero, one estimates the optical band gap Eg.

--------------------------------------------------
8. Educational use
--------------------------------------------------

This app is suitable for:

- Undergraduate and postgraduate labs on optical materials
- Quick screening of band gaps for series of glass or semiconductor samples
- Demonstrating the difference between direct and indirect transitions

You are welcome to fork and adapt it for your own teaching and research needs.

--------------------------------------------------
9. Licence
--------------------------------------------------

Insert your chosen licence text here, for example:

MIT Licence
