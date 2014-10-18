# print a list in tabular format
def make_table(table, header, justify = 'R', columnWidth = 0):
    # determine the column width if not provided
    if columnWidth == 0:
        for row in table:
            # find max column width
            for col in row:
                width = len(str(col))
                if width > columnWidth:
                    columnWidth = width

    outputStr = ""
    numCols = 0
    # make right justified columns
    for row in table:
        rowList = []
        for col in row:
            if justify == 'R':      # right justify
                rowList.append('| '+str(col).rjust(columnWidth))
            elif justify == 'L':    # justify left
                rowList.append('| '+str(col).ljust(columnWidth))
            elif justify == 'C':    # justify center
                rowList.append('| '+str(col).center(columnWidth))

        # join all the columns properly 
        string = ' '.join(rowList) + '|\n'
        if numCols == 0:
            numCols = len(string) - 1 # don't need newline character
        outputStr += string
    
    # label the columns of the table
    lList = []
    for label in header:
        lList.append('|'+str(label).center(columnWidth))
    labels = ' '.join(lList) + '|\n'
    
    # put a top and buttom banner on the table
    banner = '-'*numCols+'\n'
    outputStr = banner + labels + banner + outputStr + banner
    return outputStr
