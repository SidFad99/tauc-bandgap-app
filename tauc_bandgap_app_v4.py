
import streamlit as st
import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt

# ---------------------------------
# Helper functions
# ---------------------------------

def convert_thickness_to_cm(value, unit):
    factors = {
        "cm": 1.0,
        "mm": 0.1,
        "µm": 1e-4,
        "nm": 1e-7,
    }
    return value * factors[unit]


def alpha_from_absorbance(A, thickness_cm):
    # Beer–Lambert: A = log10(I0/I), alpha = 2.303 * A / d
    return 2.303 * A / thickness_cm


def alpha_from_transmittance(T, thickness_cm):
    # T is I/I0 (0–1). A = -log10(T)
    A = -np.log10(T)
    return alpha_from_absorbance(A, thickness_cm)


def kubelka_munk_FR(R):
    # F(R) = (1 - R)^2 / (2R), pseudo-absorption from diffuse reflectance
    eps = 1e-8
    R = np.clip(R, eps, 1 - eps)
    return (1.0 - R) ** 2 / (2.0 * R)


def linear_fit_with_Eg(hv, y):
    # Simple least-squares line: y = m x + c, Eg from x intercept
    coeffs = np.polyfit(hv, y, 1)
    m, c = coeffs
    Eg = -c / m

    # R^2 as a quick "how linear" check
    y_fit = m * hv + c
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan

    return m, c, Eg, r2


def r2_quality_text(r2):
    if np.isnan(r2):
        return "R² not defined (constant data in fit window)."
    if r2 >= 0.98:
        return "Excellent linearity (R² ≥ 0.98)."
    if r2 >= 0.95:
        return "Good linearity (0.95 ≤ R² < 0.98)."
    return "Poor linearity (R² < 0.95) – adjust the energy window."


def draw_band_diagram(Eg, title, kind="direct"):
    """
    Simple E–k band diagram.

    kind = "direct"   -> VB max and CB min at same k (vertical transition)
    kind = "indirect" -> CB min shifted in k (diagonal transition)
    """
    if not np.isfinite(Eg) or Eg <= 0:
        return None

    # k-axis (arbitrary units)
    k = np.linspace(-2.0, 2.0, 400)
    a = 0.3  # curvature

    if kind == "direct":
        # Valence band maximum at k = 0, energy = 0
        Ev = -a * k**2
        Ec = Eg + a * k**2

        k_v = 0.0
        Ev_top = 0.0
        Ec_min = Eg

        fig, ax = plt.subplots(figsize=(3, 4))
        ax.plot(k, Ev, lw=2, label="Valence band")
        ax.plot(k, Ec, lw=2, label="Conduction band")

        # Vertical transition
        ax.annotate(
            "",
            xy=(k_v, Ec_min),
            xytext=(k_v, Ev_top),
            arrowprops=dict(arrowstyle="<->", linewidth=1.5),
        )
        ax.text(
            k_v + 0.1,
            Eg / 2.0,
            f"E$_g$ ≈ {Eg:.3f} eV",
            va="center",
            fontsize=9,
        )

        ax.set_title(title + " (direct)", fontsize=11)

    else:  # indirect
        # Valence band maximum at k = 0
        Ev = -a * k**2
        # Conduction band minimum shifted in k (e.g. around k = 1)
        kc = 1.0
        Ec = Eg + a * (k - kc)**2

        k_v = 0.0        # VB max at k = 0
        k_c = kc         # CB min at k = kc
        Ev_top = 0.0
        Ec_min = Eg

        fig, ax = plt.subplots(figsize=(3, 4))
        ax.plot(k, Ev, lw=2, label="Valence band")
        ax.plot(k, Ec, lw=2, label="Conduction band")

        # Diagonal transition (photon + phonon)
        ax.annotate(
            "",
            xy=(k_c, Ec_min),
            xytext=(k_v, Ev_top),
            arrowprops=dict(arrowstyle="<->", linewidth=1.5),
        )
        ax.text(
            (k_v + k_c) / 2.0 + 0.1,
            Eg / 2.0,
            f"E$_g$ ≈ {Eg:.3f} eV",
            va="center",
            fontsize=9,
        )

        # Horizontal phonon arrow along VB (change in k)
        ax.annotate(
            "",
            xy=(k_c, Ev_top - 0.05 * Eg),
            xytext=(k_v, Ev_top - 0.05 * Eg),
            arrowprops=dict(arrowstyle="<-", linewidth=1.0),
        )
        ax.text(
            (k_v + k_c) / 2.0,
            Ev_top - 0.12 * Eg,
            "phonon (Δk)",
            ha="center",
            fontsize=8,
        )

        ax.set_title(title + " (indirect)", fontsize=11)

    ax.set_xlabel("k (wavevector, arb. units)")
    ax.set_ylabel("Energy (eV)")
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-0.4 * Eg, 1.6 * Eg)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig


# ---------------------------------
# Streamlit UI
# ---------------------------------

st.set_page_config(page_title="UV–Vis–NIR Band Gap (Tauc) App", layout="wide")

st.title("UV–Vis–NIR → Absorption Coefficient → Tauc Band Gap")
st.caption("initiative by Siddiq Fadhil, UKM Malaysia")

st.markdown(
    """
This app converts UV–Vis–NIR data to absorption-like quantities and constructs **Tauc plots**
for direct and indirect transitions. It is intended as both a research helper and a teaching tool.
"""
)

with st.expander("Show method and equations"):
    st.markdown("### Key relations")
    st.markdown("**Photon energy**")
    st.latex(r"h\nu\,[\text{eV}] = \frac{1240}{\lambda\,[\text{nm}]}")
    st.markdown("**From absorbance A (base-10) to absorption coefficient α**")
    st.latex(r"\alpha = \frac{2.303\,A}{d}")
    st.markdown("**From transmittance T = I/I_0**")
    st.latex(r"A = -\log_{10} T,\qquad \alpha = \frac{2.303\,(-\log_{10} T)}{d}")
    st.markdown("**Kubelka–Munk pseudo-absorption for diffuse reflectance**")
    st.latex(r"F(R_\infty) = \frac{(1 - R_\infty)^2}{2R_\infty}")
    st.markdown("**Tauc relation**")
    st.latex(r"(\alpha h\nu)^n \propto (h\nu - E_g)")
    st.markdown(
        r"""with \(n = 2\) for allowed direct transitions and \(n = \tfrac{1}{2}\)
        for allowed indirect transitions. The band gap \(E_g\) is obtained by fitting
        a straight line in the quasi-linear part of the Tauc plot and extrapolating
        to \((\alpha h\nu)^n = 0\)."""
    )

st.sidebar.header("Step 0 – Measurement settings")

alpha_source = st.sidebar.selectbox(
    "Quantity supplied by the spectrometer",
    [
        "Absorbance A",
        "Transmittance T",
        "Reflectance R (Kubelka–Munk)",
    ],
)

# Thickness only makes sense for A or T routes
if alpha_source in ["Absorbance A", "Transmittance T"]:
    thickness_value = st.sidebar.number_input(
        "Sample thickness", min_value=1e-6, value=0.1, format="%.6f"
    )
    thickness_unit = st.sidebar.selectbox(
        "Thickness unit", ["cm", "mm", "µm", "nm"]
    )
else:
    thickness_value = None
    thickness_unit = None

# Units for T or R
if alpha_source == "Transmittance T":
    data_unit = st.sidebar.selectbox(
        "Transmittance unit in file",
        ["fraction (0–1)", "percent (0–100)"],
    )
elif alpha_source == "Reflectance R (Kubelka–Munk)":
    data_unit = st.sidebar.selectbox(
        "Reflectance unit in file",
        ["fraction (0–1)", "percent (0–100)"],
    )
else:
    data_unit = None

st.subheader("Step 1 – Upload data")

uploaded_file = st.file_uploader(
    "Upload UV–Vis–NIR data file (CSV or TXT). The file should contain a wavelength column (nm) and one data column.",
    type=["csv", "txt"],
)

if uploaded_file is None:
    st.stop()

# Try reading with comma first, then whitespace as fallback
raw = uploaded_file.read()
try:
    df = pd.read_csv(io.BytesIO(raw))
except Exception:
    df = pd.read_csv(io.BytesIO(raw), sep=r"\s+")

st.subheader("Step 2 – Check the data preview")
st.dataframe(df.head())

# ---------------------------------
# Column selection
# ---------------------------------

st.subheader("Step 3 – Select columns")

col_wave = st.selectbox("Wavelength column (nm)", df.columns)

if alpha_source == "Absorbance A":
    col_data = st.selectbox("Absorbance A column (base-10)", df.columns)
elif alpha_source == "Transmittance T":
    col_data = st.selectbox("Transmittance T column", df.columns)
else:
    col_data = st.selectbox("Reflectance R column", df.columns)

st.caption("Tip: If the preview looks wrong, go back and check the file or column choices.")

# ---------------------------------
# Core calculations
# ---------------------------------

# Force numeric arrays for wavelength and raw signal
wavelength_nm = pd.to_numeric(df[col_wave], errors="coerce").to_numpy()
y_raw = pd.to_numeric(df[col_data], errors="coerce").to_numpy()

# Sanity checks for raw units
if alpha_source == "Absorbance A":
    if (y_raw < -0.01).any() or (y_raw > 6).any():
        st.warning("Absorbance values are outside the typical range (0–6). Check units or baseline.")
elif alpha_source == "Transmittance T":
    if data_unit == "fraction (0–1)":
        if (y_raw < 0).any() or (y_raw > 1).any():
            st.warning("Some transmittance values are outside 0–1. Did you mean percent?")
    else:
        if (y_raw < 0).any() or (y_raw > 100).any():
            st.warning("Some transmittance values are outside 0–100%. Check your data.")
elif alpha_source == "Reflectance R (Kubelka–Munk)":
    if data_unit == "fraction (0–1)":
        if (y_raw < 0).any() or (y_raw > 1).any():
            st.warning("Some reflectance values are outside 0–1. Did you mean percent?")
    else:
        if (y_raw < 0).any() or (y_raw > 100).any():
            st.warning("Some reflectance values are outside 0–100%. Check your data.")

# Photon energy (not yet sorted)
hv_eV = 1240.0 / wavelength_nm  # photon energy in eV (hc/e ≈ 1240 eV·nm)

# Compute alpha-like quantity
if alpha_source == "Absorbance A":
    d_cm = convert_thickness_to_cm(thickness_value, thickness_unit)
    alpha_like = alpha_from_absorbance(y_raw, d_cm)

elif alpha_source == "Transmittance T":
    if data_unit == "percent (0–100)":
        T = y_raw / 100.0
    else:
        T = y_raw
    d_cm = convert_thickness_to_cm(thickness_value, thickness_unit)
    alpha_like = alpha_from_transmittance(T, d_cm)

else:  # Reflectance with Kubelka–Munk
    if data_unit == "percent (0–100)":
        R = y_raw / 100.0
    else:
        R = y_raw
    alpha_like = kubelka_munk_FR(R)  # pseudo-absorption

# Data cleaning report
nan_count = np.isnan(alpha_like).sum()
neg_count = np.sum(alpha_like < 0)

alpha_like = np.nan_to_num(alpha_like, nan=0.0)
alpha_like[alpha_like < 0] = 0.0

st.subheader("Step 4 – Calculated quantities")

df["hv_eV"] = hv_eV
df["alpha_like"] = alpha_like  # cm^-1 for A/T route; F(R) for KM route
df["alpha_hv"] = df["alpha_like"] * df["hv_eV"]
df["tauc_direct"] = df["alpha_hv"] ** 2
df["tauc_indirect"] = np.sqrt(np.clip(df["alpha_hv"], a_min=0, a_max=None))

# Sort by energy for nicer plots
df.sort_values("hv_eV", inplace=True)
df.reset_index(drop=True, inplace=True)

st.dataframe(
    df[["hv_eV", "alpha_like", "alpha_hv", "tauc_direct", "tauc_indirect"]].head(10)
)

st.markdown(
    f"""
**Data summary**  
- Number of points: {len(df)}  
- Wavelength range: {np.nanmin(wavelength_nm):.1f}–{np.nanmax(wavelength_nm):.1f} nm  
- Energy range: {df["hv_eV"].min():.2f}–{df["hv_eV"].max():.2f} eV  
- NaN values replaced by zero: {nan_count}  
- Negative α values clipped to zero: {neg_count}
"""
)

# ---------------------------------
# Plots in tabs
# ---------------------------------

st.subheader("Step 5 – Inspect spectra and Tauc plots")

tab_raw, tab_direct, tab_indirect = st.tabs(
    ["Raw spectrum", "Tauc (direct)", "Tauc (indirect)"]
)

with tab_raw:
    st.markdown("**Raw spectrum**")
    try:
        raw_df = pd.DataFrame({"Wavelength_nm": wavelength_nm, "Signal": y_raw})
        raw_df = raw_df.dropna().sort_values("Wavelength_nm")
        if len(raw_df) > 0:
            st.line_chart(raw_df.set_index("Wavelength_nm"))
        else:
            st.write("No valid numeric data for raw spectrum.")
    except Exception:
        st.write("Raw spectrum plot not available – check that the chosen columns are numeric.")

with tab_direct:
    st.markdown("**Direct allowed transition:** y = (αhν)² vs hν")
    st.line_chart(df.set_index("hv_eV")[["tauc_direct"]])

with tab_indirect:
    st.markdown("**Indirect allowed transition:** y = (αhν)¹ᐟ² vs hν")
    st.line_chart(df.set_index("hv_eV")[["tauc_indirect"]])

# ---------------------------------
# Linear regression and Eg
# ---------------------------------

st.subheader("Step 6 – Choose energy windows and extract band gaps")

e_min = float(df["hv_eV"].min())
e_max = float(df["hv_eV"].max())

st.markdown(
    "Slide the handles to cover the **straight part** of each Tauc plot near the absorption edge."
)

Eg_d = None
Eg_i = None
r2_d = np.nan
r2_i = np.nan
fit_ranges = {}

col_d, col_i = st.columns(2)

with col_d:
    st.markdown("### Direct band gap fit (y = (αhν)²)")
    e1_d, e2_d = st.slider(
        "Energy window for direct fit (eV)",
        min_value=e_min,
        max_value=e_max,
        value=(max(e_min, e_min + 0.2), max(e_min + 0.4, e_max - 0.2)),
        step=0.01,
        key="direct_window",
    )
    mask_d = (df["hv_eV"] >= e1_d) & (df["hv_eV"] <= e2_d)
    n_points_d = int(mask_d.sum())
    st.caption(f"{n_points_d} points used in this direct fit.")

    if n_points_d > 2:
        x_d = df.loc[mask_d, "hv_eV"].to_numpy()
        y_d = df.loc[mask_d, "tauc_direct"].to_numpy()
        m_d, c_d, Eg_d, r2_d = linear_fit_with_Eg(x_d, y_d)
        fit_ranges["direct"] = (e1_d, e2_d)

        st.write(f"Estimated direct band gap: **Eg ≈ {Eg_d:.3f} eV**")
        st.text(f"Fit line: y = {m_d:.3e}·x + {c_d:.3e}")
        st.text(f"R² (direct) = {r2_d:.3f}")
        st.caption(r2_quality_text(r2_d))

        # Detailed plot with fit line and intercept
        fig_d, ax_d = plt.subplots()
        ax_d.plot(df["hv_eV"], df["tauc_direct"], label="Tauc data")
        ax_d.plot(x_d, m_d * x_d + c_d, label="Linear fit")
        if np.isfinite(Eg_d):
            ax_d.axvline(Eg_d, linestyle="--", label="Eg (x-intercept)")
        ax_d.set_xlabel("hν (eV)")
        ax_d.set_ylabel("(αhν)²")
        ax_d.set_title("Direct Tauc plot with linear fit")
        ax_d.legend()
        st.pyplot(fig_d)
    else:
        st.info("Please choose a wider range for the direct fit (at least 3 points).")

with col_i:
    st.markdown("### Indirect band gap fit (y = (αhν)¹ᐟ²)")
    e1_i, e2_i = st.slider(
        "Energy window for indirect fit (eV)",
        min_value=e_min,
        max_value=e_max,
        value=(max(e_min, e_min + 0.2), max(e_min + 0.4, e_max - 0.2)),
        step=0.01,
        key="indirect_window",
    )
    mask_i = (df["hv_eV"] >= e1_i) & (df["hv_eV"] <= e2_i)
    n_points_i = int(mask_i.sum())
    st.caption(f"{n_points_i} points used in this indirect fit.")

    if n_points_i > 2:
        x_i = df.loc[mask_i, "hv_eV"].to_numpy()
        y_i = df.loc[mask_i, "tauc_indirect"].to_numpy()
        m_i, c_i, Eg_i, r2_i = linear_fit_with_Eg(x_i, y_i)
        fit_ranges["indirect"] = (e1_i, e2_i)

        st.write(f"Estimated indirect band gap: **Eg ≈ {Eg_i:.3f} eV**")
        st.text(f"Fit line: y = {m_i:.3e}·x + {c_i:.3e}")
        st.text(f"R² (indirect) = {r2_i:.3f}")
        st.caption(r2_quality_text(r2_i))

        fig_i, ax_i = plt.subplots()
        ax_i.plot(df["hv_eV"], df["tauc_indirect"], label="Tauc data")
        ax_i.plot(x_i, m_i * x_i + c_i, label="Linear fit")
        if np.isfinite(Eg_i):
            ax_i.axvline(Eg_i, linestyle="--", label="Eg (x-intercept)")
        ax_i.set_xlabel("hν (eV)")
        ax_i.set_ylabel("(αhν)¹ᐟ²")
        ax_i.set_title("Indirect Tauc plot with linear fit")
        ax_i.legend()
        st.pyplot(fig_i)
    else:
        st.info("Please choose a wider range for the indirect fit (at least 3 points).")

# ---------------------------------
# Summary metrics, band diagram, and methods text
# ---------------------------------

st.subheader("Step 7 – Summary, band diagram, and report text")

metric_cols = st.columns(2)
with metric_cols[0]:
    if Eg_d is not None and np.isfinite(Eg_d):
        st.metric("Direct Eg (allowed)", f"{Eg_d:.3f} eV")
    else:
        st.metric("Direct Eg (allowed)", "N/A")
with metric_cols[1]:
    if Eg_i is not None and np.isfinite(Eg_i):
        st.metric("Indirect Eg (allowed)", f"{Eg_i:.3f} eV")
    else:
        st.metric("Indirect Eg (allowed)", "N/A")

st.markdown("### Band-gap visualisation (direct vs indirect in E–k space)")
st.caption(
    "Simplified band diagrams. For a direct gap, the valence-band maximum and conduction-band minimum "
    "are at the same k, so the transition is vertical. For an indirect gap, the conduction-band minimum "
    "is shifted in k and the transition involves both a change in energy and a change in momentum (phonon)."
)

bd_cols = st.columns(2)
if Eg_d is not None and np.isfinite(Eg_d):
    with bd_cols[0]:
        fig_bd_d = draw_band_diagram(Eg_d, "Direct band gap", kind="direct")
        if fig_bd_d is not None:
            st.pyplot(fig_bd_d)
if Eg_i is not None and np.isfinite(Eg_i):
    with bd_cols[1]:
        fig_bd_i = draw_band_diagram(Eg_i, "Indirect band gap", kind="indirect")
        if fig_bd_i is not None:
            st.pyplot(fig_bd_i)

# Auto-generated methods paragraph
meas_desc = {
    "Absorbance A": "transmission-mode UV–Vis–NIR absorbance (A)",
    "Transmittance T": "transmission-mode UV–Vis–NIR transmittance (T)",
    "Reflectance R (Kubelka–Munk)": "diffuse reflectance UV–Vis–NIR data converted using the Kubelka–Munk remission function F(R)",
}[alpha_source]

if alpha_source in ["Absorbance A", "Transmittance T"]:
    thickness_text = f"The sample thickness was {thickness_value:g} {thickness_unit}."
else:
    thickness_text = "The sample was measured in diffuse reflectance mode (optically thick layer)."

direct_range_text = (
    f"{fit_ranges['direct'][0]:.2f}–{fit_ranges['direct'][1]:.2f} eV"
    if "direct" in fit_ranges
    else "a selected region near the edge"
)
indirect_range_text = (
    f"{fit_ranges['indirect'][0]:.2f}–{fit_ranges['indirect'][1]:.2f} eV"
    if "indirect" in fit_ranges
    else "a selected region near the edge"
)

Eg_direct_text = (
    f"{Eg_d:.3f} eV (R² = {r2_d:.3f})" if Eg_d is not None and np.isfinite(Eg_d) else "not determined"
)
Eg_indirect_text = (
    f"{Eg_i:.3f} eV (R² = {r2_i:.3f})" if Eg_i is not None and np.isfinite(Eg_i) else "not determined"
)

methods_text = f"""
UV–Vis–NIR optical measurements were carried out and the data were exported as {meas_desc}.
{thickness_text}

Wavelength λ (nm) was converted to photon energy hν (eV) using hν = 1240/λ.
For transmission data, the absorption coefficient α was obtained from the base-10 absorbance A
using α = 2.303A/d, where d is the sample thickness in cm. For diffuse reflectance data,
the Kubelka–Munk remission function F(R) = (1 − R)²/(2R) was used as a pseudo-absorption term.

Tauc plots of (αhν)² (allowed direct transitions) and (αhν)¹ᐟ² (allowed indirect transitions)
were constructed as a function of hν. Straight lines were fitted to the quasi-linear parts of
these plots in the ranges {direct_range_text} for the direct case and {indirect_range_text} for the indirect case.
Extrapolation of the fitted lines to (αhν)ⁿ = 0 yielded optical band gaps of
E_g(direct) ≈ {Eg_direct_text} and E_g(indirect) ≈ {Eg_indirect_text}.
"""

st.markdown("**Suggested methods paragraph (you can copy–paste into a report):**")
st.text_area("Methods text", value=methods_text, height=260)

# ---------------------------------
# Download processed data
# ---------------------------------

st.subheader("Step 8 – Download processed data")

csv_out = df.to_csv(index=False)
st.download_button(
    "Download CSV with absorption and Tauc columns",
    csv_out,
    file_name="tauc_bandgap_processed.csv",
    mime="text/csv",
)

st.markdown(
    """
You can now use the downloaded file for further analysis (e.g. plotting in another program)
or keep it as a record of the processing performed here.
"""
)
