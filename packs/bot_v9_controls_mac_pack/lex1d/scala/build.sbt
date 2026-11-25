ThisBuild / organization := "io.github.nsram"
ThisBuild / version := "0.1.0-SNAPSHOT"
ThisBuild / scalaVersion := "3.3.4"
ThisBuild / licenses := Seq("MIT" -> url("https://opensource.org/licenses/MIT"))
ThisBuild / homepage := Some(url("https://github.com/nsram/lex1d"))

lazy val root = (project in file("."))
  .settings(
    name := "lex1d",
    description := "LEX: fast least-squares 1-D partitioning (reference implementation)",
    scalacOptions ++= Seq("-deprecation", "-feature", "-Wunused:all"),
    libraryDependencies += "org.scalameta" %% "munit" % "1.0.0" % Test,
    Test / fork := true,
    Test / javaOptions += "-Xmx2g",
    Test / baseDirectory := (ThisBuild / baseDirectory).value
  )
