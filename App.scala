import org.apache.spark.sql.SparkSession
import io.projectglow.Glow
import org.apache.spark.sql.types.LongType
import scala.collection.JavaConversions._
import java.io.File
import org.apache.spark.sql.functions.{ col, lit, when, concat }
import org.apache.spark.sql.Row
import org.apache.spark.sql.functions.{ explode, array_contains, array }
import scala.io.Source
import sys.process._
import org.apache.log4j.Logger
import scala.collection.mutable.Queue
import org.apache.spark.sql.functions._

object App {

  /**
   *
   */
  def main(args: Array[String]): Unit = {
    val spark = SparkSession
      .builder()
      .appName("GEL")
      .config("spark.hadoop.io.compression.codecs", "io.projectglow.sql.util.BGZFCodec")
      .getOrCreate()
    spark.conf.set("spark.driver.memory", "300g")
    val log = Logger.getLogger("this.getClass")
    import spark.implicits._
    Glow.register(spark)
    spark.sql("SET spark.sql.files.maxPartitionBytes=64217728")
    log.warn("Loading markers")
    val df2 = spark.read.format("csv").option("delimiter", ":").load("800k_to_extract_indexed.txt")
    val df3 = df2.select(df2("_c0"), df2("_c1").cast(LongType).as("_c1"), df2("_c2").cast(LongType))
    val df4 = df3.withColumn("_c1", col("_c1") - 1)

    //mongoexport mongodb://172.19.4.31 --db=gel-staging --collection=genomarkers --type=csv --fields=Location,Reference,Alternative  > genomarkers.csv
    //awk -F',' '{if($1<23 && (length($2) ==1) && (length($3)==1)) print $1":"$2$2"\n"$1":"$3$3"\n"$1":"$2$3; else if($1>22 && (length($2) ==1) && (length($3)==1)) print $1":"$2$2"\n"$1":"$3$3"\n"$1":"$2$3"\n"$1":"$2"\n"$1":"$3  }' genomarkers.csv |awk '{print $1":"NR}' | sed 's/^/chr/g' | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g'  > genomarkers_indexed.csv
    //
    val df2_al = spark.read.format("csv").option("delimiter", ":").load("800k_to_extract_indexed_alleles_gt.txt")
    val df3_al = df2_al.select(df2_al("_c0"), df2_al("_c1").cast(LongType).as("_c1"), df2_al("_c2"), df2_al("_c3").cast(LongType))
    val df4_al = df3_al.withColumn("_c1", col("_c1") - 1)

    val df_ids = spark.read.format("csv").option("delimiter", "\t").option("header", "true").load("id_conversion")
    
    val ref_genome_path="GRCh38_full_analysis_set_plus_decoy_hla.fa"

    val myListFile = Source.fromFile(args(0))
    val inputDir = args(1)
    val outputDir = args(2)
    val typean = args(3)

    var q = new Queue[String]
    val thread = new Thread {
      override def run {
        log.warn("Starting thread aws")
        q synchronized {

          for (line <- myListFile.getLines()) {
            var commandcp = "aws s3 cp " + line + " " + inputDir
            commandcp !

            q += line
            q.notify
            if (q.size > 50) {
              log.warn("Thread aws sleeping")
              Thread.sleep(60000)

            }
          }
        }

      }
    }
    thread.start()

    var qcm = new Queue[String]
    var qrm = new Queue[String]
    val threadcm = new Thread {
      override def run {
        log.warn("Thread cm starting")

        qcm synchronized {
          if (qcm.size < 1) {
            log.warn("Thread cm waiting")
            qcm.wait
          }
        }
        //  while(qcm.size<1){
        //   log.warn("Thread cm sleeping")
        //  Thread.sleep(120000)
        // }

        var myFileParts = qcm.dequeue()
        while (myFileParts != "END") {
          log.warn("Thread cm fileparts " + myFileParts)
          log.warn("Thread cm qcmsize " + qcm.size)

          val myFileName = myFileParts.replace(".vcf.gz", "").replace(".gvcf.gz", "")

          log.warn("Splitting by chr")
          var command = "for i in {1..22}; do grep -E \"chr$i[^0-9_]\" " + outputDir + "/" + myFileName + ".txt | cut -f 2 > " + outputDir + "/" + myFileName + "_chr$i;done"
          Seq("/bin/bash", "-c", command).!!
          command = "grep chrX " + outputDir + "/" + myFileName + ".txt | cut -f 2 > " + outputDir + "/" + myFileName + "_chrX"
          Seq("/bin/bash", "-c", command).!!
          command = "grep chrY " + outputDir + "/" + myFileName + ".txt | cut -f 2 > " + outputDir + "/" + myFileName + "_chrY"
          Seq("/bin/bash", "-c", command).!!

          command = "for i in {1..22}; do grep -E \"chr$i[^0-9_]\" " + outputDir + "/alleles_" + myFileName + ".txt | cut -f 2 > " + outputDir + "/alleles_" + myFileName + "_chr$i;done"
          Seq("/bin/bash", "-c", command).!!
          command = "grep chrX " + outputDir + "/alleles_" + myFileName + ".txt | cut -f 2 > " + outputDir + "/alleles_" + myFileName + "_chrX"
          Seq("/bin/bash", "-c", command).!!
          command = "grep chrY " + outputDir + "/alleles_" + myFileName + ".txt | cut -f 2 > " + outputDir + "/alleles_" + myFileName + "_chrY"
          Seq("/bin/bash", "-c", command).!!

          log.warn("Transversing file")
          command = "for i in {1..22}; do tr -d '\n'  < " + outputDir + "/" + myFileName + "_chr$i > " + outputDir + "/" + myFileName + "_tr_chr$i;done"
          Seq("/bin/bash", "-c", command).!!
          command = "tr -d '\n'  < " + outputDir + "/" + myFileName + "_chrX > " + outputDir + "/" + myFileName + "_chrX_tr"
          Seq("/bin/bash", "-c", command).!!
          command = "tr -d '\n'  < " + outputDir + "/" + myFileName + "_chrY > " + outputDir + "/" + myFileName + "_chrY_tr"
          Seq("/bin/bash", "-c", command).!!

          command = "for i in {1..22}; do tr -d '\n'  < " + outputDir + "/alleles_" + myFileName + "_chr$i > " + outputDir + "/alleles_" + myFileName + "_tr_chr$i;done"
          Seq("/bin/bash", "-c", command).!!
          command = "tr -d '\n'  < " + outputDir + "/alleles_" + myFileName + "_chrX > " + outputDir + "/alleles_" + myFileName + "_chrX_tr"
          Seq("/bin/bash", "-c", command).!!
          command = "tr -d '\n'  < " + outputDir + "/alleles_" + myFileName + "_chrY > " + outputDir + "/alleles_" + myFileName + "_chrY_tr"
          Seq("/bin/bash", "-c", command).!!

          
          command = "for i in {1..22}; do paste " + outputDir + "/finalIDST_" + myFileName + ".txt " + outputDir + "/" + myFileName + "_tr_chr$i >> " + outputDir + "/extracted_chr_$i;done"
          Seq("/bin/bash", "-c", command).!!
          command = "paste " + outputDir + "/finalIDST_" + myFileName + ".txt " + outputDir + "/" + myFileName + "_chrX_tr >> " + outputDir + "/extracted_chr_X"
          Seq("/bin/bash", "-c", command).!!
          command = "paste " + outputDir + "/finalIDST_" + myFileName + ".txt " + outputDir + "/" + myFileName + "_chrY_tr >> " + outputDir + "/extracted_chr_Y"
          Seq("/bin/bash", "-c", command).!!

          command = "for i in {1..22}; do paste " + outputDir + "/finalIDST_" + myFileName + ".txt " + outputDir + "/alleles_" + myFileName + "_tr_chr$i >> " + outputDir + "/extracted_alleles_chr_$i;done"
          Seq("/bin/bash", "-c", command).!!
          command = "paste " + outputDir + "/finalIDST_" + myFileName + ".txt " + outputDir + "/alleles_" + myFileName + "_chrX_tr >> " + outputDir + "/extracted_alleles_chr_X"
          Seq("/bin/bash", "-c", command).!!
          command = "paste " + outputDir + "/finalIDST_" + myFileName + ".txt " + outputDir + "/alleles_" + myFileName + "_chrY_tr >> " + outputDir + "/extracted_alleles_chr_Y"
          Seq("/bin/bash", "-c", command).!!

          qrm synchronized {
            qrm += myFileParts
            log.warn("Thread cm notify Thread rm")
            qrm.notify
          }
          qcm synchronized {
            if (qcm.size < 1) {
              log.warn("Thread cm waiting")
              qcm.wait
            }
          }
          log.warn("Dequeue de qcm")
          myFileParts = qcm.dequeue()

        }
        qrm synchronized {
          log.warn("Thread cm put END in Thread rm")
          qrm += "END"
          qrm.notify
        }

      }
    }

    threadcm.start()
    val threadrm = new Thread {
      override def run {

        log.warn("Thread rm starting")
        qrm synchronized {
          if (qrm.size < 1) {
            log.warn("Thread rm waiting")
            qrm.wait
          }
        }

        

        var myFileParts = qrm.dequeue()
        while (myFileParts != "END") {
          log.warn("Thread rm fileparts " + myFileParts)
          log.warn("Thread rm qcmsize " + qrm.size)

          val myFileName = myFileParts.replace(".vcf.gz", "").replace(".gvcf.gz", "")
          log.warn("Removing file:" + inputDir + "/" + myFileParts)
          val removecp = "rm " + inputDir + "/" + myFileParts
          Seq("/bin/bash", "-c", removecp).!!

          log.warn("Removing directory:" + outputDir + "/" + myFileName)
          var removedr = "rm -r " + outputDir + "/" + myFileName + "*"
          Seq("/bin/bash", "-c", removedr).!!
          removedr = "rm -r " + outputDir + "/alleles_" + myFileName + "*"
          Seq("/bin/bash", "-c", removedr).!!
          
          qrm synchronized {
            if (qrm.size < 1) {
              log.warn("Thread rm waiting")
              qrm.wait
            }
          }
          log.warn("Dequeue de qcm")
          myFileParts = qrm.dequeue()

        }

      }
    }

    threadrm.start()

    Thread.sleep(20000)

    while (q.size > 0) {
      val line = q.dequeue
      log.warn("Queue size is " + q.size)
      log.warn("Starting analysis of:" + line)
      val myFileParts = line.split("/").last
      val myFileName = myFileParts.replace(".vcf.gz", "").replace(".gvcf.gz", "")
      //Extract the ID from the file name
      val myFNAL = myFileName.split("\\.")
      var myTempId = ""
      myFNAL.take(1).foreach(row => myTempId = row.toString)
      log.warn(myTempId)
      val file = new File(inputDir + "/" + myFileParts)
      log.warn("Loading file in Glow:" + file.getAbsolutePath())
      //Read gVCF file
      val df = spark.read.format("vcf").load(file.getAbsolutePath())
      //val df = spark.read.format("vcf").load("/mnt/vol1/LP3000067-DNA_H01.genome.vcf.gz")
      if (typean == "ind-pos" || typean == "ind-all") {
        //Filter reference positions
        //Multiallelic
        //Filters PASS
        //ReferenceAlleles N
        //Genotypes not empty
        //Select cols contigName,start
        //Join with markers
        
        //Replace with 0s and 1s
        //Normalize
        //Extract reference positions and explode by END  -- NOT
        //val refDF = df.where(size(col("alternateAlleles")) === 0).withColumn("newstart", sequence(col("start"), col("end") - 1)).withColumn("newstart2", explode($"newstart")).drop("newstart").drop("start").withColumnRenamed("newstart2", "start").select("contigName", "start", "end", "referenceAllele", "alternateAlleles", "filters", "genotypes")
        //
        //val unionDF = maDF.drop("INFO_OLD_MULTIALLELIC").union(refDF)
        //First the filters because multiallelics not output this field
        var joinDF = df.where(array_contains($"filters", "PASS"))
         
        joinDF = Glow.transform("split_multiallelics", joinDF).select("contigName", "start", "end", "referenceAllele", "alternateAlleles", "genotypes")
        joinDF=joinDF.where("referenceAllele!='N'")
          .withColumn("first_genotype", explode($"genotypes")).filter(array_contains($"first_genotype.calls", 0) or array_contains($"first_genotype.calls", 1) or array_contains($"first_genotype.calls", 2))
          joinDF.cache() 
        //if (typean == "ind-pos") {
          //Right join in position with markers file
          var joinDFpos = joinDF.select(col("contigName"), col("start"), col("end"),col("referenceAllele"),col("alternateAlleles"),col("contigName").alias("jd"))
            .join(df4.withColumnRenamed("_c0", "contigName").withColumnRenamed("_c1", "start").hint("broadcast"), Seq("contigName", "start"), "right")
            .withColumn("jdb", when(col("jd").isNull, "0").otherwise("1")).orderBy(asc("_c2"))
          log.warn("Writing file in directory:" + outputDir + "/" + myFileName)
          //Comment because is costly and is only interesting for indels
          //joinDFpos=Glow.transform("normalize_variants", joinDFpos, Map("reference_genome_path" -> ref_genome_path))
          joinDFpos.select("contigName", "jdb", "_c2").distinct().orderBy(asc("_c2")).repartition(1).write.format("csv").option("delimiter", "\t").option("header", "false").save(outputDir + "/" + myFileName)

        //} else if (typean == "ind-all") {
          var joinDFall = joinDF.filter(length(col("referenceAllele")) < 2)
            .where(size(col("alternateAlleles")) === 0 || length(col("alternateAlleles")(0)) < 2)
            .withColumn("genotype", col("first_genotype.calls"))

            .withColumn("gt", when(array_contains(col("genotype"), 0) && !array_contains(col("genotype"), 1), concat(col("referenceAllele"), col("referenceAllele")))
              .otherwise(when(array_contains(col("genotype"), 1) && !array_contains(col("genotype"), 0), concat(col("alternateAlleles")(0), col("alternateAlleles")(0)))
                .otherwise(when(array_contains(col("genotype"), 0) && array_contains(col("genotype"), 1), concat(col("referenceAllele"), col("alternateAlleles")(0)))
                  .otherwise(when((col("contigName") === "chrX" || col("contigName") === "chrY") && array_contains(col("genotype"), 0) && size(col("genotype")) < 2, col("referenceAllele"))
                    .otherwise(when((col("contigName") === "chrX" || col("contigName") === "chrY") && array_contains(col("genotype"), 1) && size(col("genotype")) < 2, col("alternateAlleles")(0)))))))

            //.withColumn("output",when(col("referenceAllele")===col("ALLELES") && array_contains(col("genotype"),0),"1")
            //    .otherwise(when(col("referenceAllele")=!=col("ALLELES")&&array_contains(col("genotype"),1),"1")
            //       .otherwise(when(col("referenceAllele")===col("ALLELES")&& !(array_contains(col("genotype"),0)),"0"))))

            .select(col("contigName"), col("start"), col("end"),col("referenceAllele"),col("alternateAlleles"),col("gt"), col("contigName").alias("jd"))
            .join(df4_al.withColumnRenamed("_c0", "contigName").withColumnRenamed("_c1", "start").withColumnRenamed("_c2", "gt"), Seq("contigName", "start", "gt"), "right")
            //.withColumn("jdb",when(col("jd").isNull,"0").otherwise("1")).orderBy(asc("_c2"))
            .withColumn("jd", when(col("jd").isNull, "0").otherwise("1"))
            
            //Comment because is costly and is only interesting for indels
            //joinDFall=Glow.transform("normalize_variants", joinDFall, Map("reference_genome_path" -> ref_genome_path))
          log.warn("Writing file in directory:" + outputDir + "/alleles_" + myFileName)
          joinDFall.select("contigName", "jd", "_c3").distinct().orderBy(asc("_c3")).repartition(1).write.format("csv").option("delimiter", "\t").option("header", "false").save(outputDir + "/alleles_" + myFileName)
        //}

        /* Old version with the explode of reference positions
         *
        val joinDF2 = unionDF
          .where(array_contains($"filters", "PASS"))
          .withColumn("first_genotype", explode($"genotypes"))
          .filter(array_contains($"first_genotype.calls", 0) or array_contains($"first_genotype.calls", 1) or array_contains($"first_genotype.calls", 2))
          .drop("first_genotype")
          .select(col("contigName"), col("start"), col("contigName").alias("jd"), col("end"))
          .join(df4, col("_c0").equalTo(col("contigName")).and(expr("array_contains(sequence(start),end)),_c1")), "right")
          //.join(df4,  col("_c0").equalTo(col("contigName")).and(col("_c1").geq(col("start"))).and(col("_c1").leq(col("end")))   , "right")
          .withColumn("jdb", when(col("jd").isNull, "0").otherwise("1")).orderBy(asc("_c2"))
        //.withColumnRenamed("_c0", "contigName").withColumnRenamed("_c1", "start")
         *
         */

      } else if (typean == "agg-pos") {

      } else if (typean == "agg-all") {

      } else {

      }

      log.warn("Copying file")
      val myDirectory = new File(outputDir + "/" + myFileName)
      for (file <- myDirectory.listFiles if file.getName endsWith "csv") {
        "cp " + file + " " + outputDir + "/" + myFileName + ".txt" !
      }
      val myDirectoryall = new File(outputDir + "/alleles_" + myFileName)
      for (file <- myDirectoryall.listFiles if file.getName endsWith "csv") {
        "cp " + file + " " + outputDir + "/alleles_" + myFileName + ".txt" !
      }
      
      

      val finalId = df_ids.where("platekey == '" + myTempId + "'").select("participant_id")
      var finalIDST = ""
      finalId.collect.foreach(row => finalIDST = row(0).toString)

      val commandIDST = "echo " + finalIDST + " >" + outputDir + "/finalIDST_" + myFileName + ".txt"
      Seq("/bin/bash", "-c", commandIDST).!!
      log.warn(finalIDST)

      qcm synchronized {
        qcm += myFileParts
        log.warn("Thread cm notified")
        qcm.notify
      }

    }
    qcm synchronized {
      log.warn("Put END in Thread cm")
      qcm += "END"
      qcm.notify
    }

  }

}