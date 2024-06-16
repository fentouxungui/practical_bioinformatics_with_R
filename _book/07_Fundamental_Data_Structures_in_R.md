--- 
always_allow_html: true
---



# Fundamental Data Structures in R

## Named Vectors

A named vector allows you to assign names to individual elements within the vector. This may seem like a small feature, but it can greatly enhance your ability to organize, manipulate, and analyze data effectively.

### What is a Named Vector?
In R, a vector can be named, meaning that each element within the vector can have a descriptive name associated with it. Think of it as a way to label your data. You can use the `names()` function to create a named vector by assigning names to the elements within the vector.

Let's start with a simple example:


```r
# Create a numeric vector
expression_vec <- c(10, 25, 30, 12, 20)

# Assign names to the vector elements
names(expression_vec) <- c("gene1", "gene2", "gene3", "gene4", "gene5")

# View the named vector
expression_vec
```

```
## gene1 gene2 gene3 gene4 gene5 
##    10    25    30    12    20
```

As you can see, each element now has a corresponding name, making it easier to identify and work with specific values in the vector.

### Using Names to Subset a Vector
One of the key advantages of named vectors is the ability to subset them based on the names. This allows you to access specific elements of the vector easily.

For example, if you want to select only the values associated with "gene1" and "gene3," you can do so like this:


```r
# Subset the vector using names
selected_genes <- expression_vec[c("gene1", "gene3")]

# View the selected genes
selected_genes
```

```
## gene1 gene3 
##    10    30
```

### Named Vectors as Dictionary-Like Structures

Unlike some programming languages like Python, R doesn't have a built-in dictionary data structure. However, named vectors can serve as a similar tool because they establish a one-to-one mapping between names and values.

Imagine you have gene expression measurements across various conditions that you want to analyze. You can create a named vector like this:


```r
# Create a named vector from scratch
expression <- c(normal = 2.3, treated = 5.1, resistant = 3.8, sensitive = 1.2)

# View the named vector
expression
```

```
##    normal   treated resistant sensitive 
##       2.3       5.1       3.8       1.2
```

Now, if you want to retrieve the expression value for the "resistant" condition, you can do so by using the name as follows:


```r
# Access the expression value for "resistant" condition
expression["resistant"]
```

```
## resistant 
##       3.8
```

### Real-life example
In single-cell RNAseq analysis, you have tens of thousands of cells in the experiment and you cluster the cells into different clusters (cell types). Usually, you get the cluster id as numeric numbers:1,2,3,4,5,...

After you annotate the clusters with marker genes, you can give the clusters a specific name.


```r
cells<- c("cell1", "cell2", "cell3", "cell4", "cell5", "cell6")

# the clustering algorithum assign each cell to cluster 1-4 
clusters<- c(1, 2, 3, 2, 1, 4)

# we have annotated the clusters 1-4 as follows:
annotation<- c("CD4T", "CD8T", "NK", "B cell")

# create a named vector
names(annotation)<- c(1,2,3,4)

annotation
```

```
##        1        2        3        4 
##   "CD4T"   "CD8T"     "NK" "B cell"
```

We can use this named vector to re-annotate the original cells


```r
annotation[clusters]
```

```
##        1        2        3        2        1        4 
##   "CD4T"   "CD8T"     "NK"   "CD8T"   "CD4T" "B cell"
```
We can then combine the cells and the new annotation as a data frame (which we will cover later).


```r
data.frame(cell= cells, cell_annotation= annotation[clusters])
```

```
##    cell cell_annotation
## 1 cell1            CD4T
## 2 cell2            CD8T
## 3 cell3              NK
## 4 cell4            CD8T
## 5 cell5            CD4T
## 6 cell6          B cell
```

You see that the named vector can be very useful as a dictionary.

### Reordering and Subsetting with Names
One of the remarkable features of named vectors is their flexibility. Even if the order of names differs from another vector, you can still use the names to reorder or subset the vector effectively.

Let's say you have another vector representing fold changes in gene expression:


```r
# Create a vector for fold change of gene expression
folds <- c(resistant = 1.1, sensitive = 4.2, normal = 1.3, treated = 2.1)

# View the fold change vector
folds
```

```
## resistant sensitive    normal   treated 
##       1.1       4.2       1.3       2.1
```

Notice that the order of conditions is different from the expression vector. However, you can use the names from the expression vector to reorder or subset the fold change vector:


```r
# Reorder the fold change vector using names from the expression vector
folds[names(expression)]
```

```
##    normal   treated resistant sensitive 
##       1.3       2.1       1.1       4.2
```

### Conclusion
Named vectors may seem simple, but they offer a valuable way to organize and manipulate your data, serving as a powerful tool for tasks like subsetting, indexing, and organizing information. This beginner-friendly guide has introduced you to the concept of named vectors and demonstrated their practical use in real-world scenarios. As you delve deeper into R, you'll find that mastering this fundamental feature will greatly enhance your data analysis capabilities.

## Lists
Lists in R are versatile data structures that can hold various types of elements, including both numeric values and character strings, and even elements of different lengths. In this section, we will explore the concept of lists, understand how to create and manipulate them, and learn different ways to access their elements.

### Introduction to Lists
Unlike matrices and vectors, which can only store elements of the same type, lists can accommodate a mix of numeric and character elements with different lengths. This flexibility makes lists a powerful tool for organizing and storing complex data structures.

### Creating a List
Imagine you are conducting a series of biological experiments over multiple days. For each day, you collect various data, including numeric measurements (e.g., gene expression levels) and descriptive summaries (e.g., experimental conditions). Lists can be used to store this data efficiently. Each list element can represent a day's data, containing both numeric vectors and character strings to capture the diverse information for each experiment:


```r
results <- list(p_values = c(0.01, 0.56),
                summary = "Non-significant",
                experiments = c("day1", "day2"))
results
```

```
## $p_values
## [1] 0.01 0.56
## 
## $summary
## [1] "Non-significant"
## 
## $experiments
## [1] "day1" "day2"
```

When we print the `results` list, you will notice that each element is preceded by a `$` sign, indicating its name.

### Accessing List Elements
You can access list elements using various methods:

#### Using $ Notation
To access an element by its name, you can use the `$` notation:


```r
results$experiments
```

```
## [1] "day1" "day2"
```

```r
results$summary
```

```
## [1] "Non-significant"
```

#### Using Single Brackets
Using single brackets `[ ]` will return a list:


```r
results[1]
```

```
## $p_values
## [1] 0.01 0.56
```

####Slicing the List
You can slice a list by specifying a range of indices:


```r
results[1:2]
```

```
## $p_values
## [1] 0.01 0.56
## 
## $summary
## [1] "Non-significant"
```

#### Using `[[ ]]` Double Brackets

To access the actual element rather than a list, use double brackets `[[ ]]`:


```r
results[[1]]
```

```
## [1] 0.01 0.56
```

```r
results[["p_values"]]
```

```
## [1] 0.01 0.56
```

#### Using Names
You can also access elements by their names using single brackets and double brackets:


```r
# returns a list
results["p_values"]
```

```
## $p_values
## [1] 0.01 0.56
```

```r
# return the element 
results[["p_values"]]
```

```
## [1] 0.01 0.56
```

**Tip**: In RStudio, when you type `$` after the list, it will display a dropdown menu showing all the available names of the list elements, making it easier to access them.

### Conclusions

Lists are essential for managing diverse data structures efficiently in R. They allow you to store and organize data with different characteristics, making them a valuable tool for data manipulation and analysis in various applications.

## Dataframes

Data frames are one of the fundamental data structures in R, and they play a pivotal role in various data analysis tasks. In the realm of biology research, data frames are exceptionally useful for organizing and manipulating biological data efficiently. In this tutorial, we will explore the basics of data frames, including their creation, column access, and subsetting.

### What is a Data Frame?
A data frame is a two-dimensional tabular data structure in R. It is similar to a spreadsheet or a database table, where rows and columns intersect to form a grid-like structure. In a data frame, each row can hold different types of data, such as numbers, text, or factors, but all columns must have the same number of rows and each column should have the same data type.

### Creating a Data Frame

To create a data frame, you can use the `data.frame()` function. Here's an example of creating a simple data frame for patient data in a biology study.


```r
patient_data <- data.frame(
  id = c("P1", "P2"),
  age = c(42, 36),
  status = c("normal", "tumor")
)
```

In this example, we have three vectors (id, age, and status) that are combined into a data frame. The variable names of the vectors become the column names of the data frame. You can see the resulting data frame by simply typing `patient_data`.


```r
patient_data
```

```
##   id age status
## 1 P1  42 normal
## 2 P2  36  tumor
```

### Accessing Columns
You can access individual columns within a data frame using the `$` sign, which is similar to accessing elements in a list. For instance:


```r
patient_data$id
```

```
## [1] "P1" "P2"
```

This command will retrieve the `id` column from the `patient_data` data frame.

Similarly, you can access the `status` column by using `patient_data$status`.

### Subsetting a Data Frame
Subsetting a data frame means extracting specific rows and columns based on your analysis needs. You can use square brackets `[row_indices, column_indices]` to subset or slice the data frame, similar to subsetting a matrix.

For example, to select the first row and the first and third columns from the patient_data data frame:


```r
patient_data[1, c(1, 3)]
```

```
##   id status
## 1 P1 normal
```

This command returns a subset of the data frame with the `id` and `status` columns for the first row.

### Conclusion
In biology research, data frames are invaluable for organizing and analyzing diverse datasets. Imagine you have a dataset with patient information, where each row represents a patient and columns contain attributes like age, gender, and diagnosis. You could use data frames to filter patients by specific criteria, calculate summary statistics, or create visualizations to gain insights into your biological data.

Data frames in R serve as the backbone for handling and exploring data in a structured and meaningful way, making them an essential tool for any biologist or data scientist working in the field of life sciences.

## How to Represent Categorical Data in R? Understanding Factors

In R, factors are a fundamental data type used to represent categorical data. Factors are akin to biological categories, such as gene types, experimental conditions, or mutation classes. They provide a structured way to handle and analyze qualitative data. In this tutorial, we will explore factors in R, starting from their creation, manipulation, and practical applications.

### Creating Factors
Let's begin with a hypothetical scenario involving phenotypes: wild type (WT), mutant (Mut), and heterozygote (Het). We can create a factor called `phenotypes` to represent these categories:


```r
phenotypes <- factor(c("WT", "Mut", "Mut", "WT", "Het"))
```

When you print `phenotypes`, you'll see the categories along with their associated levels:


```r
phenotypes
```

```
## [1] WT  Mut Mut WT  Het
## Levels: Het Mut WT
```

In this example, the factor phenotypes has three levels: "Het," "Mut," and "WT."

### Changing levels
By default, the levels are ordered alphabetically. You can use the relevel to change it


```r
relevel(phenotypes, ref = "WT")
```

```
## [1] WT  Mut Mut WT  Het
## Levels: WT Het Mut
```

Now, the `WT` becomes the base level.

### Converting Factors

Factors can be converted to characters or numeric values. When you convert them to numbers, R assigns a numerical value to each category based on their order in the factor levels. The lowest level gets assigned the number 1, and the highest level gets the highest number.

#### As character

This R command converts the phenotypes factor into a character vector. In other words, it changes the factor's categorical labels into text.


```r
as.character(phenotypes)
```

```
## [1] "WT"  "Mut" "Mut" "WT"  "Het"
```

#### Converting to numbers

This R command converts the `phenotypes` factor into a numeric vector. It assigns a numerical value to each category based on their order in the factor levels.


```r
as.numeric(phenotypes)
```

```
## [1] 3 2 2 3 1
```

In this case, "WT" corresponds to 3, "Mut" to 2, and "Het" to 1.

### Customizing Factor Levels
You can customize the order of factor levels to control how they are converted to numbers. Suppose you want "Mut" to be assigned the lowest value. You can specify the desired order when creating the factor:


```r
phenotypes <- factor(c("WT", "Mut", "Mut", "WT", "Het"), 
                     levels= c("Mut", "Het", "WT"))
```

Now, when you convert it to numbers, "Mut" will be 1:


```r
as.numeric(phenotypes)
```

```
## [1] 3 1 1 3 2
```

### Creating Factors from Character Vectors
You can also create a factor from a character vector and explicitly specify the desired factor levels. For example:


```r
phenotypes <- c("WT", "Mut", "Mut", "WT", "Het")
phenotypes <- factor(phenotypes, levels= c("Mut", "Het", "WT"))

phenotypes
```

```
## [1] WT  Mut Mut WT  Het
## Levels: Mut Het WT
```

This approach allows you to control the factor levels directly.

### Subsetting Factors
When you subset a factor, the factor levels are preserved. For instance:


```r
phenotypes[1:2]
```

```
## [1] WT  Mut
## Levels: Mut Het WT
```

Even though we only extracted elements 1 and 2, the factor levels "Mut," "Het," and "WT" are still present.

### Removing Unused Levels
Sometimes, you may want to remove unused factor levels to keep your data clean. You can achieve this using the `droplevels()` function. For example:


```r
droplevels(phenotypes[1:2])
```

```
## [1] WT  Mut
## Levels: Mut WT
```
In this case, the "Het" level was removed because it was not present in the subset.

### Conclusion
Factors are invaluable in various fields of research, including biology and genetics, where they are used to represent and analyze categorical data like gene types, mutation classes, or experimental conditions. Understanding factors and their manipulation is essential for accurate data analysis and interpretation.

Statistical modeling uses factors heavily. For example, when doing differential gene expression analysis, the factor levels determine which group is used as the baseline. Also, when plotting figures, the order of the geometric objects (e.g, the boxes in the boxplot) is ordered by the factor levels. Changing the factor levels will change the look of the plots. We will see an example in the later data visulization lecture.

In summary, factors in R provide a structured way to work with categorical data, allowing you to control the order of levels and efficiently manage and analyze information. Mastery of factors is a crucial skill for any data scientist or researcher working with categorical data in R.

## Section Complete 
Congratulations on completing the section!

You've learned about named vectors, lists, data frames, and factors. Understanding these basics is essential for working with data effectively. This knowledge will be valuable as you tackle more advanced data analysis tasks.

As you move on to the next section, remember to use our Q&A section and comments for any questions you have. Let's keep exploring and learning together in our R programming journey!


                     
