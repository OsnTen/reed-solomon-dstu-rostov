import random
import sys
from random import randint
s=0
sox=0
koll = 0
match len( sys.argv):
    case 2:
        s= sys.argv[1]
        koll=int(input("Введите количество ошибок: "))
        sox=input("Введите путь для сохранения: ")
        
    case 3:
        s=sys.argv[1]
        koll=int(sys.argv[2])
        sox=input("Введите путь для сохранения: ")
    case 4:
        s=sys.argv[1]
        koll=int(sys.argv[2])
        sox=sys.argv[3]
    case _:
        s=input("Введите путь для открытия файла: ")
        koll=int(input("Введите количество ошибок: "))
        sox=input("Введите путь для сохранения: ")
f = open(s, 'rb')


tMesg = f.read()


errospisok =[]
class ReedSolomon:
# поля Галуа

    __GFEXP = [0] * 65536
    __GFLOG = [0] * 256

    def __init__(self):  # инициализация конструктора
        self.__GFEXP[0] = 1
        byteValu = 1
        for bytePos in range(1 ,255):
            byteValu <<= 1
            if (byteValu & 0x100):
                byteValu ^= 0x11d

            self.__GFEXP[bytePos] = byteValu
            self.__GFLOG[byteValu] = bytePos

        for bytePos in range(255 ,65536):
            self.__GFEXP[bytePos] = self.__GFEXP[bytePos - 255]

    def __gfMult(self, argX, argY):
        if ((argX == 0) or (argY == 0)):
            byteValu = 0
        else:
            byteValu = self.__GFLOG[argX]
            byteValu += self.__GFLOG[argY]
            byteValu = self.__GFEXP[byteValu]
        return (byteValu)

    def __gfDivi(self, argX, argY):
        if (argY == 0):
            raise ZeroDivisionError()
        if (argX == 0):
            byteValu = 0
        else:
            byteValu = self.__GFLOG[argX] - self.__GFLOG[argY]
            byteValu += 255
            byteValu = self.__GFEXP[byteValu]
        return (byteValu)

    def _gfPolyAdd(self, polyA, polyB):
        polySum = [0] * max(len(polyA), len(polyB))

        for polyPos in range(0, len(polyA)):
            polySum[polyPos + len(polySum) - len(polyA)] = polyA[polyPos]

        for polyPos in range(0, len(polyB)):
            polySum[polyPos + len(polySum) - len(polyB)] ^= polyB[polyPos]
        return (polySum)

    def _gfPolyMult(self, polyA, polyB): # перемножение полиномов
        polyProd = len(polyA) + len(polyB) - 1
        polyProd = [0] * polyProd
        for posB in range(0, len(polyB)):
            for posA in range(0, len(polyA)):
                polyProd[posA + posB] ^= self.__gfMult(polyA[posA], polyB[posB])
        return (polyProd)

    def _gfPolyScale(self, argPoly, argX):  # полиноминальное масштабирование
        polyVal = [0] * len(argPoly)
        for polyPos in range(0, len(argPoly)):
            polyVal[polyPos] = self.__gfMult(argPoly[polyPos], argX)
        return (polyVal)

    def _gfPolyEval(self, argPoly, argX): # полиномиальная оценка
        byteValu = argPoly[0]

        for polyPos in range(1, len(argPoly)):
            tempValu = self.__gfMult(byteValu, argX)
            tempValu = tempValu ^ argPoly[polyPos]
            byteValu = tempValu
        return (byteValu)

    def _rsGenPoly(self, errSize): # представление полинома
        polyValu = [1]

        for polyPos in range(0, errSize):
            tempVal = [1, self.__GFEXP[polyPos]]
            polyValu = self._gfPolyMult(polyValu, tempVal)
        return (polyValu)

    def RSEncode(self, argMesg, errSize):
        polyGen = self._rsGenPoly(errSize)

        outBuffer = (len(argMesg) + errSize)
        outBuffer = [0] * outBuffer

        for mesgPos in range(0, len(argMesg)):
            mesgChar = argMesg[mesgPos]
            outBuffer[mesgPos] = (mesgChar)

        for mesgPos in range(0, len(argMesg)): # кодирование
            mesgChar = outBuffer[mesgPos]
            if (mesgChar != 0):
                for polyPos in range(0, len(polyGen)):
                    tempValu = self.__gfMult(polyGen[polyPos], mesgChar)
                    outBuffer[mesgPos + polyPos] ^= tempValu
        for mesgPos in range(0, len(argMesg)):
            mesgChar = argMesg[mesgPos]
            outBuffer[mesgPos] = (mesgChar)
        return (outBuffer)

    def _rsSyndPoly(self, argCode, errSize): # Декодирование, вычисление синдромного многочлена
        polyValu = [0] * errSize
        for errPos in range(0, errSize):
            byteValu = self.__GFEXP[errPos]
            polyValu[errPos] = self._gfPolyEval(argCode, byteValu)
        return (polyValu)

    def _rsForney(self, polySynd, eraseLoci, errSize):  # алгоритм Форни для вычисления значения ошибки
        polyValu = list(polySynd)

        for posI in range(0, len(eraseLoci)):
            termX = errSize - 1 - eraseLoci[posI]
            termX = self.__GFEXP[termX]
            for posJ in range(0, len(polyValu) - 1):
                termY = self.__gfMult(polyValu[posJ], termX)
                termY ^= polyValu[posJ + 1]
                polyValu[posJ] = termY
            polyValu.pop()
        return (polyValu)

    def _rsFindErr(self, errLoci, errSize): # поиск местонахождения ошибки
        errPoly = [1]
        tempPoly = [1]

        for posSynd in range(0, len(errLoci)): # алгоритм Берлекмпа-Месси
            tempPoly.append(0)
            termSynd = errLoci[posSynd]

            for posErr in range(1, len(errPoly)):
                termPoly = errPoly[len(errPoly) - posErr - 1]
                termPoly = self.__gfMult(termPoly, errLoci[posSynd - posErr])
                termSynd ^= termPoly

            if (termSynd != 0):
                if (len(tempPoly) > len(errPoly)):
                    tNewP = self._gfPolyScale(tempPoly, termSynd)
                    tempPoly = self._gfPolyScale(errPoly, self.__gfDivi(1, termSynd))
                    errPoly = tNewP
                tempValu = self._gfPolyScale(tempPoly, termSynd)
                errPoly = self._gfPolyAdd(errPoly, tempValu)
        errCount = len(errPoly) - 1 # подсчет количество ошибок
        if ((errCount * 2) > len(errLoci)):
            print("Слишком много ошибок")
            return (None)
        else:
            print("Количество ошибок: ", errCount)
        errList = [] # вычисление нулей полинома(многочлен ошибок E)
        for errPos in range(0, errSize):
            errZed = self._gfPolyEval(errPoly, self.__GFEXP[255 - errPos])
            if (errZed == 0):
                errZed = errSize - errPos - 1
                errList.append(errZed)
        global errorspisok
        errList = errospisok
        if (len(errList) != errCount):
            print("Не удалось найти ошибки")
            print(errList)
            return (None)
        else:
            return (errList)
    def _rsCorrect(self, argCode, polySynd, errList): # исправляем ошибки и стирания
        polyLoci = [1]
        for errPos in range(0, len(errList)):
            errTerm = len(argCode) - errList[errPos] - 1
            errTerm = self.__GFEXP[errTerm]
            polyLoci = self._gfPolyMult(polyLoci, [errTerm, 1])
        errEval = polySynd[0:len(errList)] # оценка многочлена локаторов ошибки
        errEval.reverse()
        errEval = self._gfPolyMult(errEval, polyLoci)
        tMark = len(errEval) - len(errList)
        errEval = errEval[tMark:len(errEval)]
        errLoci = polyLoci[len(polyLoci) % 1 : len(polyLoci) : 2] # многочлен локаторов ошибок отнимаем четные многочлены(корни)
        for errPos in range(0, len(errList)): # корректировка
            errByte = errList[errPos] - len(argCode) + 256
            errByte = self.__GFEXP[errByte]
            errValu = self._gfPolyEval(errEval, errByte)
            errAdj = self.__gfMult(errByte, errByte)
            errAdj = self._gfPolyEval(errLoci, errAdj)
            mesgByte = self.__gfMult(errByte, errAdj)
            mesgByte = self.__gfDivi(errValu, mesgByte)
            argCode[errList[errPos]] ^= mesgByte
        return (argCode)
    def RSDecode(self, argCode, errSize): # основная часть декодирования
        codeBuffer = list(argCode)
        eraseCount = []
        for codePos in range(0, len(codeBuffer)):
            if (codeBuffer[codePos] < 0):
                codeBuffer[codePos] = 0
                eraseCount.append(codePos)
        if (len(eraseCount) > errSize):
            print("Много стираний")
            return (None)
        polySynd = self._rsSyndPoly(codeBuffer, errSize) # поиск синдромного многочлена
        if (max(polySynd) == 0):
            print("В сообщении нет ошибок")
            return (codeBuffer)
        errLoci = self._rsForney(polySynd, eraseCount, len(codeBuffer))
        errList = self._rsFindErr(errLoci, len(codeBuffer)) # ищем ошибки
        if (errList == None):
            print("Не удалось найти никаких ошибок")
            return (None)
        else:
            print("Обнаруженные ошибки ", errList)
        outMesg = self._rsCorrect(codeBuffer, polySynd, (eraseCount + errList)) # исправляем ошибки
        return (outMesg)
def baggent(tCode,tSize,k):
    global errospisok
    if k%2 !=0:
        for i in range(k):
            tekch = randint(0,tSize)
            errospisok.append(tekch)
            tCode[tekch] = randint(0,255)
    else:
        for i in range(k-1):
            tekch = randint(0,tSize)  
            errospisok.append(tekch)
            tCode[tekch] = randint(0,255)
    return tCode
def bit_fin(frames,tMesg,tSize):
    global errospisok
    if len(errospisok) > tSize/2:
        for i in range(tSize):
            tBit[i] = (tMesg[i])
        frames = bytes(tBit)
        return frames
    else:
        return frames

fooRS = ReedSolomon()
tSize = len(tMesg)
frames = tMesg
f.close()    
tCode = fooRS.RSEncode(tMesg, tSize)
tBit=[0] * tSize
for i in range(tSize):
    tBit[i] = str(tCode[i])
kod = open("koding.txt", "wt")
kodd = ' '.join(tBit)
kod.write(kodd)
kod.close()
tCode = baggent(tCode,tSize,koll)
kods = open("bag.txt", "wt")
for i in range(tSize):
    tBit[i] = str(tCode[i])
bagkod = ' '.join(tBit)
kods.write(bagkod)
kods.close()
tMesg = fooRS.RSDecode(tCode, tSize)
frames = bit_fin(frames,tMesg,tSize)
out_file = open(sox, "wb")
lol =out_file.write(frames)
out_file.close()
print("Количесвто обработанных байт: ", lol)
