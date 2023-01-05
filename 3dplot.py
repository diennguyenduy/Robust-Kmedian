import plotly.express as px
import pandas as pd

df = pd.read_csv('Kmedians_result.csv')

fig = px.scatter_3d(df,
                    x='Age',
                    y='Annual Income (k$)',
                    z='Spending Score (1-100)',
                    color='Cluster',
                    symbol='Cluster')

fig.show()