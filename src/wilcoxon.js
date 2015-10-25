var quicksort = require('./sort/quicksort')
var heapsort = require('./sort/heapsort')
var mergesort = require('./sort/mergesort')
var insertionsort = require('./sort/insertionsort')

function normalZ(z) {

    var x = z
    var b = [0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429]
    var p = 0.2316419
    var t = 1/(1+p*x)
    var fact = t
    var sum = 0
    for(var i=0; i <= b.length - 1; i++) {
        sum += b[i]*fact
        fact *= t
    }
    p = 2 * sum * Math.exp(-x*x/2.0) / (Math.sqrt(2*Math.PI))
    return p
}

var exp = module.exports

exp.wilcoxon = function(a1, a2, a1a2) {
    
    a1 = a1.slice(0)
    heapsort(a1)
    if (a2) {
        a2 = a2.slice(0)
        quicksort(a2)
    }
    
    if (!a1a2) {
        a1a2 = []
        // TODO sorting not necessary, lists already sorted
        for (var x=0; x<a1.length; x++) a1a2.push(a1[x])
        for (var y=0; y<a2.length; y++) a1a2.push(a2[y])
        insertionsort(a1a2)
    } else {
        a1a2 = a1a2.slice(0)
        quicksort(a1a2)
    }
    
    var nA = a1.length
    var n = a1a2.length
    var nB = n - nA
    var maxSum = n * (n + 1) / 2
    var h0 = maxSum / 2

    var previous = -Number.MAX_VALUE
    var start = 0
    var totalRank = []
    for (var i = 0; i<a1a2.length; i++) {
        if (a1a2[i]==previous) {
            var meanRank = (start + i + 2) / 2
            for (var j = start; j<=i; j++) {
                totalRank[j] = meanRank
            }
        } else {
            totalRank[i] = i + 1
            previous = a1a2[i]
            start = i
        }
    }

    var shortest = a1
    if (a2 && a2.length < a1.length) shortest = a2
    var nShortest = shortest.length
    var w = 0

    var r1 = 0
    var index = 0
    for (var a = 0; a < shortest.length; a++) {
        for (var j = index; j < a1a2.length; j++) {
            if (shortest[a] === a1a2[j]) {
                r1 += (j+1)
                break
            }
            index++
        }
        w += totalRank[index]
    }

    var nZ = nShortest
    if (w > h0) {
        nZ = n - nShortest
        w = maxSum - w
    }

    var uA = r1 - nA * (nA + 1) / 2
    var auc = uA / (nA * nB)

    var k = nA
    var n0 = n
    var permutations = 1
    while (k > 0) {
        k--
        n0--
        permutations *= n0/k
    }
    
    if (permutations >= 25000 || shortest.length >= 10) {
        var continuity = 0.5
        if (w>=h0) continuity = -0.5
        var z = Math.abs((w + continuity - nZ * (n + 1) / 2) / Math.sqrt(nA * nB * (n + 1) / 12))
        var p = Math.max(Number.MIN_VALUE, normalZ(z))

        return {p: p, auc: auc}
    }

    return {p: -1, auc: -1}
}

exp.wilcoxonRanked = function(a1, a2) {

    if (a2.length < a1.length) {
        var temp = a1
        a1 = a2
        a2 = temp
    }

    var nA = a1.length
    var nB = a2.length
    var n = nA + nB
    var maxSum = n * (n + 1) / 2
    var h0 = maxSum / 2
    var w = 0

    var r1 = 0
    var index = 0
    for (var a = 0; a < nA; a++) {
        r1 += a1[a]
        // TODO w on total ranks in case of ties
        w += a1[a]
    }

    var nZ = nA
    if (w>h0) {
        nZ = n - nA
        w = maxSum - w
    }

    var uA = r1 - nA * (nA + 1) / 2
    var auc = uA / (nA * nB)

    var k = nA
    var n0 = n
    var permutations = 1
    while (k > 0) {
        k--
        n0--
        permutations *= n0/k
    }
    if (permutations >= 25000 || a1.length >= 10) {
        var continuity = 0.5
        if (w>=h0) continuity = -0.5
        var z = Math.abs((w + continuity - nZ * (n + 1) / 2) / Math.sqrt(nA * nB * (n + 1) / 12))
        var p = Math.max(Number.MIN_VALUE, normalZ(z))

        return {p: p, auc: auc}
    }

    return {p: -1, auc: -1}
}
