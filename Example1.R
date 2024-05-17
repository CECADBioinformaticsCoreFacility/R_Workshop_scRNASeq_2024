# Load data
data("iris")

# I'm not sure what's in the iris data. I'm printing it below to take a look
head(iris)


## Data analysis
# ----------------------------------------------------------------------------

# Below, we calculate the average Sepal.Length across all Species in the mtcars data frame.
print(paste0("Average of Sepal Length : ", mean(iris$Sepal.Length)))

# Here, we also plot Petal.Length against Sepal.Length.
png("scatter.png")
plot(iris$Petal.Length, iris$Sepal.Length, pch=19, col=1)
dev.off()

# Typical barplot to check the mean of Petal.Length across species
png("bar.png")
df = aggregate(iris[,1:4], by = list(iris$Species), FUN = mean)
barplot(Petal.Length~Group.1, data = df,xlab = c('Species'), ylab = c('Petal Length'), col=c(1:3))
dev.off()

# Here, we plot Petal.Length across Species to check the distribution
png("box.png")
boxplot(Petal.Length ~ Species, data = iris, col = c(1:3))
dev.off()

# Subset data
setosa <- iris[iris$Species == "setosa", ]
versicolor <- iris[iris$Species == "versicolor", ]

# Run t-test
t.test(x = setosa$Petal.Length, y = versicolor$Petal.Length)

