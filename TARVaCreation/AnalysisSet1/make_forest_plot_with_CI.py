import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Typography and Style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# 2. Data Preparation
df = pd.read_csv('Fisher_Exact_unique_modified_patterns.csv')
df['gene_category'] = df['gene_category'].replace('AC-enriched', 'CN-enriched')

# Apply Haldane–Anscombe correction (+0.5 to zero counts) for stability
cols = ["AD_unique", "AD_common", "CN_unique", "CN_common"]
for col in cols:
    df[col] = df[col].astype(float).replace(0, 0.5)

# Calculate Odds Ratio (OR) and 95% Confidence Interval (CI)
df["OR_adj"] = (df["AD_unique"] * df["CN_common"]) / (df["AD_common"] * df["CN_unique"])
se = np.sqrt(1/df["AD_unique"] + 1/df["AD_common"] + 1/df["CN_unique"] + 1/df["CN_common"])
df["CI_low"] = np.exp(np.log(df["OR_adj"]) - 1.96 * se)
df["CI_high"] = np.exp(np.log(df["OR_adj"]) + 1.96 * se)

# 3. Organize Sequence (Headers above sections)
df_ad = df[df['gene_category'] == 'AD-enriched'].sort_values('OR_adj', ascending=False)
df_cn = df[df['gene_category'] == 'CN-enriched'].sort_values('OR_adj', ascending=False)

plot_sequence = []
plot_sequence.append("AD_HEADER")
for _, row in df_ad.iterrows():
    plot_sequence.append(row.to_dict())

plot_sequence.append("SPACER")

plot_sequence.append("CN_HEADER")
for _, row in df_cn.iterrows():
    plot_sequence.append(row.to_dict())

# Reverse to plot top-down
plot_sequence = plot_sequence[::-1]

# 4. Multi-Panel Layout
fig = plt.figure(figsize=(14, 10))
ax_labels = fig.add_axes([0.05, 0.1, 0.18, 0.8])
ax_pvals  = fig.add_axes([0.24, 0.1, 0.10, 0.8])
ax_ci     = fig.add_axes([0.35, 0.1, 0.14, 0.8])  # Added text panel for 95% CI values
ax_plot   = fig.add_axes([0.51, 0.1, 0.44, 0.8])

for ax in [ax_labels, ax_pvals, ax_ci]:
    ax.axis('off')
    ax.set_ylim(-1, len(plot_sequence))

ax_labels.text(0, len(plot_sequence), "Variant Call-type", weight='bold', va='bottom', fontsize=13)
ax_pvals.text(0.5, len(plot_sequence), "P-value", weight='bold', va='bottom', ha='center', fontsize=13)
ax_ci.text(0.5, len(plot_sequence), "95% CI", weight='bold', va='bottom', ha='center', fontsize=13)

# 5. Populate Data Rows with Significance Bolding
color_ad = '#cc3311' # red
color_cn = '#0077bb' # blue

for i, item in enumerate(plot_sequence):
    if isinstance(item, str):
        if item == "AD_HEADER":
            ax_labels.text(0, i, "AD-enriched Genes", weight='bold', color=color_ad, va='center', fontsize=12)
            ax_plot.axhline(i, color=color_ad, alpha=0.08, linewidth=12)
        elif item == "CN_HEADER":
            ax_labels.text(0, i, "CN-enriched Genes", weight='bold', color=color_cn, va='center', fontsize=12)
            ax_plot.axhline(i, color=color_cn, alpha=0.08, linewidth=12)
    else:
        # Determine Color and Significance
        cat_color = color_ad if item['gene_category'] == 'AD-enriched' else color_cn
        p_val = item['p_value']
        
        # BOLD LOGIC: p < 0.05
        txt_weight = 'bold' if p_val < 0.05 else 'normal'
        
        # Variant Name
        ax_labels.text(0.1, i, item['mod_type'], va='center', fontsize=11, fontweight=txt_weight)
        
        # P-value String (bold if significant)
        p_str = f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.3f}"
        ax_pvals.text(0.5, i, p_str, va='center', ha='center', fontsize=11, fontweight=txt_weight)
        
        # CI Value String (bold if significant)
        ci_str = f"[{item['CI_low']:.2f}, {item['CI_high']:.2f}]"
        ax_ci.text(0.5, i, ci_str, va='center', ha='center', fontsize=11, fontweight=txt_weight)
        
        # Forest Plot Circles and Bars
        ax_plot.errorbar(item['OR_adj'], i, 
                         xerr=[[item['OR_adj'] - item['CI_low']], [item['CI_high'] - item['OR_adj']]],
                         fmt='o', color=cat_color, ecolor=cat_color, 
                         capsize=3, markersize=6, elinewidth=1.2, markeredgewidth=1.2)

# 6. Final Axis Styling
ax_plot.set_ylim(-1, len(plot_sequence))
ax_plot.axvline(1, color='black', linestyle='--', alpha=0.4, linewidth=1)
ax_plot.set_xscale('log')
ax_plot.set_xlabel("Odds Ratio (log scale, 95% CI)", fontsize=13, fontweight='bold')
ax_plot.set_yticks([])

plt.savefig('Figure4_forest_plot.png', dpi=300, bbox_inches='tight')
