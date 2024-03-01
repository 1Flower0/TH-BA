class PatternSearch:

    def badCharHeuristic(string,size):
        badChar = {'A':-1,'C':-1,'T':-1,'G':-1,'+':-1}
        for i in range(size):
            badChar[string[i]]=i
        return badChar

    def search(self,txt,pat):
        patternLength = len(pat)
        textLength = len(txt)
        indicies = []
        badChar = PatternSearch.badCharHeuristic(pat,patternLength)
        index = 0
        while(index <= textLength-patternLength):
            patternIndex = patternLength-1
            while patternIndex >=0 and pat[patternIndex]==txt[index+patternIndex]:
                patternIndex-=1
            if patternIndex<0:
                indicies.append(index)
                print(indicies.append(index))
                index+=(patternLength-badChar[txt[index+patternLength]] if index+patternLength<textLength else 1)
            else:
                index+=max(1,patternIndex-badChar[txt[index+patternIndex]])
        return indicies
    
    def __init__(self) -> None:
        pass