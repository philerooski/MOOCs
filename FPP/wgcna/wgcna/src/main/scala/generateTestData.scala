package com.example.wgcna

import scala.util.Random
import java.io._

object generateTestData {
    def randomGenes(m: Int): List[String] = {
        def iter(m: Int, genes: List[String]): List[String] = {
            if (m == 0) genes
            else {
                val newGene = ( Seq.fill(4)(97 + Random.nextInt(25)) 
                    map { _.toChar } ).mkString
                iter(m - 1, newGene :: genes)
            }
        }
        iter(m, List())
    }

    def main(args: Array[String]) = {
        val (m, n) = (args(0).toInt, args(1).toInt)
        val genes = randomGenes(m)
        val samples = ( (1 to n) map { "S" + _.toString } ).toList
        val expressionMatrix = for (i <- 1 to m) yield {
            ( for (j <- 1 to n) yield Random.nextDouble.toString ).toList
        }
        val dataset = for (l <- (genes zip expressionMatrix)) yield {
            (l._1 :: l._2).mkString("\t")
        }
        val header = ( "GeneNames" :: samples ).mkString("\t")
        val pw = new PrintWriter(new File("testData.tsv"), "UTF-8")
        pw.write(( header :: dataset ).mkString("\n"))
        pw.close
    }
}
