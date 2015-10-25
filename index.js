'use strict'

var descriptives = require('./src/descriptives')
var ttest = require('./src/ttest')
var wilcoxon = require('./src/wilcoxon')

module.exports = {

    mean: descriptives.mean,
    variance: descriptives.variance,
    stdev: descriptives.stdev,

    covariance: descriptives.covariance,
    correlation: descriptives.correlation,
    correlationStdNorm: descriptives.correlationStdNorm,

    standardNormalize: descriptives.standardNormalize,
    
    welch: ttest.welch,
    student: ttest.student,
    wilcoxon: wilcoxon.wilcoxon,
    wilcoxonRanks: wilcoxon.wilcoxonRanks
    
}
