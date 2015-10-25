var descriptives = require('./descriptives')
var probability = require('./probability')

var exp = module.exports

exp.student = function(a1, a2) {

    var n1 = a1.length
    var n2 = a2.length
    var mean1 = descriptives.mean(a1)
    var mean2 = descriptives.mean(a2)
    var var1 = descriptives.variance(a1, mean1)
    var var2 = descriptives.variance(a2, mean2)

    var sx1x2 = Math.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2));
    var t = (mean1 - mean2) / (sx1x2 * Math.sqrt(1 / n1 + 1 / n2));
    var df = n1 + n2 - 2;
    
    var p = t < 0 ? p = probability.studentTCDF(df, t) : probability.studentTCDF(df, -t)
    if (p === 0) p = Number.MIN_VALUE

    return {t: t, p: p, df: df}
}

exp.welch = function(a1, a2) {

    var n1 = a1.length
    var n2 = a2.length
    var mean1 = descriptives.mean(a1)
    var mean2 = descriptives.mean(a2)
    var var1 = descriptives.variance(a1, mean1)
    var var2 = descriptives.variance(a2, mean2)

    var t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);
    var df = ((var1 / n1 + var2 / n2) * (var1 / n1 + var2 / n2)) / (((var1 / n1) * (var1 / n1)) / (n1 - 1) + ((var2 / n2) * (var2 / n2)) / (n2 - 1));

    var p = t < 0 ? p = probability.studentTCDF(df, t) : probability.studentTCDF(df, -t)
    if (p === 0) p = Number.MIN_VALUE

    return {t: t, p: p, df: df}
}
