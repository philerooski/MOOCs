package com.example.wgcna

import scala.io.Source
import scala.collection.parallel.immutable.ParSeq

object wgcna {
    def readTsv(f: String): List[Array[String]] = {
        val expressionByLine = Source.fromFile(f).getLines.toList
        expressionByLine map { _.split("\t") }
    }

    def parseExpressionData(d: List[Array[String]]) = {
        val sampleIndex = d.head.tail
        val geneIndex = d.tail map { _(0) } 
        val expressionMatrix = d.tail map { _.tail map { _.toDouble } }
        (sampleIndex, geneIndex, expressionMatrix)
    }

    def residual(x: Array[Double]): Array[Double] = {
        val meanX = x.sum / x.length
        x map { _ - meanX }
    }

    def pearsonCorrelation(resX: Array[Double], resY: Array[Double]): Double = {
        val numer = ( { resX zip resY } map { case (x, y) => x*y } ).sum
        def squareSumSqrt(res: Array[Double]) = {
            scala.math.pow( (res map { scala.math.pow(_, 2) }).sum, 0.5 )
        }
        val denom = squareSumSqrt(resX) * squareSumSqrt(resY)
        numer / denom
    }

    def correlationMatrix(expressionMatrix: List[Array[Double]]): ParSeq[List[Double]] = {
        val eMatrixPar = expressionMatrix.par
        val residuals = eMatrixPar map { residual(_) }
        for (x1 <- residuals) yield {
            for (x2 <- residuals.toList) yield pearsonCorrelation(x1, x2)
        }
    }

    def main(args: Array[String]) {
        val expression = readTsv("testData.tsv")
        val (sampleIndex, geneIndex, expressionMatrix) = parseExpressionData(expression)
        val corMatrix = correlationMatrix(expressionMatrix)
        println(corMatrix.head)
        println("Hello World!!!")
    }
}
