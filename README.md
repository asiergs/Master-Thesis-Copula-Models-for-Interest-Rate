# Master Thesis Copula Models for Interest Rate
Full master thesis can be read in the following link [MasterThesis_CopulaModelProposalForInterestRateAndEquityShocksAggregationUnderSolvencyII.pdf](https://github.com/asiergs/Master-Thesis-Copula-Models-for-Interest-Rate/blob/main/MasterThesis_CopulaModelProposalForInterestRateAndEquityShocksAggregationUnderSolvencyII.pdf). The PowerPoint presentation (in spanish) can be found in the following link.
Insurance companies, due to their business nature, are exposed to multiple risks. These risks begin with the signing of an insurance contract, where the company commits to some uncertain future payments, which in some cases can be doubly uncertain in terms of frequency and severity.
However, claims risks are not the only risks that insurance companies must face. For example, there are additional risks due to the financial component needed to be provided in the insurance products, where guaranteed returns are sometimes promised, especially in long-term products such as lifetime annuities. To achieve these returns, companies must invest in the market by acquiring financial assets that provide the necessary profitability to meet their obligations, with assets varying from fixed income to equity values such as stocks or linked bonds.
The Solvency II regulation recently introduced and mandatory in many European countries since 2016 establishes a framework that considers all the risks to which an insurance company is exposed. This way all insurance companies can provide a more self-adapted solvency capital calculation in comparison to the calculation required by previous regulation. This way, the capital requirements criteria are also unified so all companies are considering the same risk scheme defined in the regulation.
The following illustration shows the risk map proposed by the regulation Solvency II:
<p align="center">
   <img src="https://raw.githubusercontent.com/asiergs/Master-Thesis-Copula-Models-for-Interest-Rate/main/solvencyIIrisksmap.svg" alt="2400"/>
</p>
When it comes to risk aggregation, the use of the standard formula provided by the regulation together with the risk correlations provided is the most common practice for insurance companies. However, there are other more complex alternatives that allow for a more appropriate aggregation of company risks, such as internal models that can use the same risk map or another if considered more appropriate.
Master thesis in Actuarial Science
2
The underlying assumptions under the standard model proposed by the regulation are very restrictive and mainly conditioned by the normality of the data. This can easily be proven to be true in some modules due to the central limit theorem, but in the case of the market module it is well known that financial series are not governed by normality. In these cases, the standard formula cannot properly reflect the real market behavior and missing the impact the market shocks could have in the solvency of the insurance companies.
The standard formula makes simplifications for the sake of simplicity so the results can be easily applied and interpreted by any risk manager. However, when trying to determine the real trade-off this simplification provides, there is hardly any literature, if not any, about a reliable comparison of this simplified approach and a statistically backed method.
Additionally, the recent high inflation scenarios have proven the weaknesses of the interest rate model and shocks provided by the regulation and a new model calibration that properly reflects the market behavior is required for the proper aggregation modelling.
