## 29 April 2025
On plane from Denver --> Vancouver
Dothis: Code line 50 of cleanTS_2025.R. Double check with Lizzie if I should really do the cleaning step of selection only one phenophase when the same one was recorded by the different people.

**Multiple_FirstY**: Indicates whether there are multiple series for the phenophase on the organism within the selected time period. A value of "0" indicates there is only one series, and a value of "1" indicates there are more than one series within the selected time period. When the selected time period consists of more than one year, calculation of multiple series is separate for each 12-month time span starting from the selected start date. Note that the first "yes" record for a 12-month time span during the selected time period might not be preceded by a "no" record, but the first "yes" for subsequent series within the same 12-month time span is always separated from the previous series by at least one "no" record.
**Multiple_Observers:** Indicates whether multiple observers contributed records to the series. A value of "0" indicates only one observer contributed records, and a value of "1" indicates more than one observer contributed records. Multiple observers can contribute observations at a shared group site (http://www.usanpn.org/nn/groups/shared-site).

**NumYs_in_Series**: The total number of days in the series with a "yes" record. 


I should also revise how I calculate GDD. For now, in cleanTS 2025, its calculated from budburst to leaf drop, or when absent, leaf coulouring.