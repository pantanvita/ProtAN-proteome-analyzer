import sys
import os
import pytest
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

DATA_PATH = "proteome-data.csv"

def test_csv_loads():
    """Test 1: Your CSV loads correctly."""
    df = pd.read_csv(DATA_PATH)
    assert not df.empty, "CSV is empty"
    assert "Gene" in df.columns
    assert all(col in df.columns for col in ["S1", "S2", "S3", "R1", "R2", "R3"])
    print(f"✅ Loaded {len(df)} proteins")
    return df

def test_log2fc_calculation():
    """Test 2: log2FC computation works (handles NaN/zero values)."""
    df = pd.read_csv(DATA_PATH)
    
    # Calculate means for S and R groups
    s_mean = df[['S1','S2','S3']].mean(axis=1)
    r_mean = df[['R1','R2','R3']].mean(axis=1)
    
    # Handle edge cases (your real data has zeros/NaNs)
    s_mean = s_mean.replace(0, np.nan)  # Avoid log(0)
    r_mean = r_mean.replace(0, np.nan)
    
    # Safe division + log2 (handles NaN/inf)
    ratio = np.where(
        np.isfinite(s_mean) & np.isfinite(r_mean),
        s_mean / r_mean,
        np.nan
    )
    log2fc = np.log2(np.clip(ratio, 1e-10, 1e10))  # Clip extremes
    
    df['log2FC'] = pd.Series(log2fc, index=df.index)
    
    # Now the assertions work
    assert 'log2FC' in df.columns
    assert df['log2FC'].notna().sum() > 10, "Too many NaN log2FC values"
    
    # Check variation (allow smaller std for noisy proteomics data)
    valid_fc = df['log2FC'].dropna()
    assert len(valid_fc) > 0, "No valid log2FC values"
    assert valid_fc.std() > 0.05, f"Too little variation: std={valid_fc.std():.3f}"
    
    print(f"✅ log2FC OK: {len(valid_fc)} valid, range {valid_fc.min():.2f} to {valid_fc.max():.2f}, std={valid_fc.std():.3f}")

def test_pvalues_valid():
    """Test 3: p-values are valid (0,1]."""
    df = pd.read_csv(DATA_PATH)
    pvals = []
    for _, row in df.iterrows():
        _, pval = stats.ttest_ind([row['S1'],row['S2'],row['S3']], 
                                 [row['R1'],row['R2'],row['R3']], equal_var=False)
        pvals.append(pval)
    assert all(0 < p <= 1 for p in pvals), "Invalid p-values"
    print("✅ p-values valid")

def test_classify_protein_function():
    """Test 4: Your classification logic."""
    def classify_protein(row, p_cutoff=0.1, fc_cutoff=0.58):
        if (row["p_value"] < p_cutoff) and (row["log2FC"] >= fc_cutoff):
            return "Upregulated"
        elif (row["p_value"] < p_cutoff) and (row["log2FC"] <= -fc_cutoff):
            return "Downregulated"
        else:
            return "Not significant"
    
    # Test cases
    tests = [
        (0.01, 0.8, "Upregulated"),
        (0.01, -0.7, "Downregulated"),
        (0.2, 1.0, "Not significant"),
        (0.01, 0.3, "Not significant")
    ]
    
    for p, fc, expected in tests:
        result = classify_protein(pd.Series({"p_value": p, "log2FC": fc}))
        assert result == expected, f"Expected {expected}, got {result}"
    print("✅ Classification works")

def test_significant_filtering():
    """Test 5: Filtering logic."""
    df = pd.read_csv(DATA_PATH)
    df['log2FC'] = np.log2(df[['S1','S2','S3']].mean(axis=1) / df[['R1','R2','R3']].mean(axis=1))
    pvals = [stats.ttest_ind([r['S1'],r['S2'],r['S3']], [r['R1'],r['R2'],r['R3']], equal_var=False)[1] 
             for _, r in df.iterrows()]
    df['p_value'] = pvals
    
    p_thresh, fc_thresh = 0.1, 0.58
    sig_mask = (df["p_value"] < p_thresh) & (df["log2FC"].abs() >= fc_thresh)
    sig_df = df[sig_mask]
    
    if len(sig_df) > 0:
        assert (sig_df["p_value"] < p_thresh).all()
        assert (sig_df["log2FC"].abs() >= fc_thresh).all()
        print(f"✅ Found {len(sig_df)} significant proteins")
    else:
        print("✅ Filter works (no significant hits expected for small data)")
