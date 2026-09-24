plugins {
    application
}

group = "jackSolver"
version = "0.1.0"

java {
    toolchain {
        languageVersion = JavaLanguageVersion.of(21)
    }
}

repositories {
    mavenCentral()
}

dependencies {
    // CERN Colt matrix library (pulls in concurrent:concurrent transitively)
    implementation("colt:colt:1.2.0")

    testImplementation("junit:junit:4.13.2")
}

application {
    mainClass = "jackSolver.JackRHFSolver"
}

tasks.test {
    useJUnit()
    testLogging {
        events("passed", "skipped", "failed")
        exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.FULL
    }
}
