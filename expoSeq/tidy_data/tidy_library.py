def calcFraction(df,x, target_column):
    return df[target_column] / df[target_column].sum()


def summarizeDuplicates(df, target_column, dependency):
    column = df[target_column].sum()
    df[target_column] = column
    df = df.drop_duplicates(dependency)
    return df

mapping_funcs = {
    "calcFraction": calcFraction,
    "summarizeDuplicates": summarizeDuplicates,
}
