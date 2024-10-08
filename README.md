# Prion Aggregate (Propagon) Replication Rates & Asymmetric Transmission
Code for estimation of replication rate and asymmetric distribution using experimental aggregate count 
publication: <a href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010107">PLOS Computational Biology Article</a>

Prion aggregates, cell division dynamics, and transmission of aggregates. <br /><img src="figures/Fig0.png" alt="Prion Aggregate Replication Dynamics." width="550" />

<strong>Included:</strong>
<ul>
  <li><b><em>propagon_data_raw</em></b>: Folder containing experimental aggregate count datasets for six prion variants.<br></li>
  <li><b><em>propagon_data_filtered_iqr</em></b>: Folder containing aggregate counts data filtered for outliers.<br></li>
  <li><b><em>propagon_data_down_sampled</em></b>: Folder containing aggregate counts data downsampled in time.<br></li>
  <li><b><em>Example_Simulate_Data.m</em></b>: The script shows how to generate simulated data using the functions in the <b><em>/code/</em></b> folder.<br></li>
  <li><b><em>Example_Parameter_Estimation.m</em></b>: The script shows how to apply the adaptive Metropolis algorithm to simulated and experimental data.<br></li>
</ul>
<br />

<strong>Mathematical Modeling and Replication Rate Estimation</strong>
<ol> 
   <li>Prion aggregate count experiments.<br /><img src="figures/Fig1.png" alt="Prion Aggregate Count Experiments." width="550" /> </li>
   <li>Structured model, aggregate distribution, and cell division.<br /><img src="figures/Fig2.png" alt="Diagram of Structured Model of Prion Aggregate Distribution." width="550" /> </li>
   <li>Prion aggregate counts for six prion variants.<br /><img src="figures/Fig3.png" alt="Prion Aggregate Count for Six Prion Variants." width="550" /> </li>
   <li>Adaptive Metropolis estimates of prion aggregate replication rates. <br /><img src="figures/Fig4.png" alt="Prion replication rates estimation" width="550" /> </li>
   <li>Adaptive Metropolis estimates of prion aggretate transmission bias. <br /><img src="figures/Fig5.png" alt="Prion replication rates estimation" width="550" /> </li><br />
  <em>Note: The appearance of increasing variance over time is a visual artifact from viewing iterations on a logarithmic scale.</em>
</ol>
