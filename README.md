# genstats
Statistical tests and descriptives for science

Mathematics regarding statistical tests are adapted from [Colt](https://dst.lbl.gov/ACSSoftware/colt/) and [JSci](http://jsci.sourceforge.net/) Java libraries. 

## Descriptives

* arithmetic mean
* sample variance
* standard deviation
* sample covariance
* sample correlation

## Statistical tests

* Student's t-test (equal variances t-test)
* Welch's t-test (unequal variances t-test)
* Wilcoxon test (Mann-Whitney U test)

T-tests return a p-value, a t-statistic and degrees of freedom. Wilcoxon test returns a p-value and area under the curve (AUC). All tests are one-tailed. Multiply the returned p-value by two to get a two-tailed p-value.

## Usage example
```javascript
var genstats = require('genstats')

var a1 = [], a2 = []
for (var i = 0; i < 10000; i++) {
    a1.push(Math.random())
    a2.push(Math.random() - 0.05)
}

console.log('Descriptives')

console.log('mean\t\t' + genstats.mean(a1))
console.log('variance\t' + genstats.variance(a1))
console.log('stdev\t\t' + genstats.stdev(a1))
console.log('stdev^2\t\t' + genstats.stdev(a1) * genstats.stdev(a1))

console.log('covariance\t' + genstats.covariance(a1, a2))
console.log('correlation\t' + genstats.correlation(a1, a2))

console.log()
console.log('Statistical tests')

// Student's t-test                                                                                                                                                                
console.log('Student\'s t-test', genstats.student(a1, a2))

// Welch's t-test (unequal variances)                                                                                                                                              
console.log('Welch\'s t-test', genstats.welch(a1, a2))

// Wilcoxon test (Mann-Whitney U)                                                                                                                                                  
console.log('Wilcoxon', genstats.wilcoxon(a1, a2))

```

## Test
```bash
$ npm test
```
