# Code for the KSMS proposal in 2019A

Most of BJ's work happens in notebooks/Simulator SB.ipynb. I reccommend you copy that and start messing around in your own notebook.

Adjust the subsamples defined as classes in `sim.py` to tweak the survey.

Example:

```python
## Perhaps we want to turn this sample into an M-dwarf selection
class SampleClose(Sample):
    name = 'close'
    description = '10 closest dwarfs that are bighter than V < 13'
    def in_sample(self):
        df = self.df.copy()
        b = pd.Series(False, index=self.df.index) 
        idx = self.df.query('_DEJ2000 > 0 and Vmag < 13 and Rs < 1.5').sort_values(by='Dist').iloc[:10].index
        b.loc[idx] = True
        return b
```
Change the `idx = ` line to 
```python
idx = self.df.query('_DEJ2000 > 0 and Vmag < 13 and Rs < 1.5').sort_values(by='Dist').iloc[:30].index
```
in order to include the closest 30 stars, instead of 10. The possibilities are endless!


The process for tweaking is as follows:
1. Clone the repo and create a branch for yourself.
1. Tweak the survey as you wish.
1. Submit a pull request

Have fun!
