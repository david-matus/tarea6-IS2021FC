import numpy as np
import matplotlib.pyplot as plt
import os

def CrearDirectorio(folder):
    '''
    Crea una carpeta con el nombre dado en la ubicación donde se está ejecutando el programa.
    Si existe una carpeta con el mismo nombre, crea una con la forma nombre (1), etc.
    Entradas:
    folder: str: nombre deseado de la carpeta
    Salidas:
    subfolder: str: nombre con el que fue creada la carpeta
    '''
    try: #primero se intenta guardar con un nombre sencillo
        os.mkdir(folder)
        subfolder = folder
    except:
        filepend = True
        trycount = 1
        while filepend:
            try: #se añade una cola al nombre de carpeta si la carpeta ya existe
                os.mkdir(folder+' ('+str(trycount)+')')
                filepend = False
                subfolder = folder+' ('+str(trycount)+')'
            except:
                trycount += 1
    return subfolder #el nombre de la carpeta creada

def InicializarPoblación(tamañoPoblación, nGenes):
    '''
    Crea una matriz de tamaño tamañoPoblaciónxnGenes que contiene los cromosomas iniciales.
    Cada gen es creado aleatoriamente, pero no se repiten dentro del cromosoma.
    Entradas:
    tamañoPoblación: int: cantidad de cromosomas deseados
    nGenes: int: cantidad de genes que debe tener cada cromosoma
    Salidas:
    población: nxm matrix: lista con todos los cromosomas
    '''
    población = []
    for i in range(tamañoPoblación): #cada ciclo es un nuevo cromosoma
        cromosoma = []
        genesInsertados = 0
        while genesInsertados<nGenes: #cada ciclo se intenta añadir un gen
            sinRepetir = True
            contadorGenes = 0
            nuevoGen = np.random.randint(0,nGenes)
            while (contadorGenes<genesInsertados and sinRepetir):
                sinRepetir = (cromosoma[contadorGenes]!=nuevoGen)
                contadorGenes += 1
            if sinRepetir: #sólo se añade el gen si no está repetido
                cromosoma.append(nuevoGen)
                genesInsertados += 1
        población.append(cromosoma)
    return población

def EvaluarPoblación(población,ciudadesXY):
    '''
    Calcula el valor de evaluación de cada uno de los cromosomas en la población. La evaluación
    es el inverso de la distancia total recorrida por la secuencia del cromosoma.
    Entradas:
    población: nxm matrix: matriz que contiene los cromosomas
    ciudadesXY: 2xn matrix: matriz [X,Y] que contiene las coordenadas de cada ciudad
    Salidas:
    puntuaciones: 1xn numpy array: lista de puntuaciones de cada cromosoma
    '''
    Totalcromosomas = len(población)
    Totalgenes = len(población[0])
    puntuaciones = np.zeros(Totalcromosomas)
    for icromosoma in range(Totalcromosomas):
        puntuacionActual = 0
        cromosomaMod = np.copy(población[icromosoma])
        cromosomaActual = np.append(cromosomaMod,cromosomaMod[0])
        for igen in range(Totalgenes):
            a = cromosomaActual[igen]
            b = cromosomaActual[igen+1]
            x1, y1 = ciudadesXY[0][a], ciudadesXY[1][a]
            x2, y2 = ciudadesXY[0][b], ciudadesXY[1][b]
            puntuacionActual += np.sqrt((x1-x2)**2+(y1-y2)**2)
        puntuaciones[icromosoma] = 1/puntuacionActual
    return puntuaciones

def OperadorMutación(individuo):
    '''
    Modifica un cromosoma invirtiendo la posición de dos genes aleatoriamente.
    Entradas:
    individuo: 1xn matrix: cromosoma a mutar
    Salidas:
    cromosomaMutado: 1xn matrix: cromosoma con genes invertidos
    '''
    gen1 = np.random.randint(0,len(individuo))
    busqueda = True
    while busqueda:
        gen2 = np.random.randint(0,len(individuo))
        busqueda = gen1==gen2
    cromosomaMutado = np.copy(individuo)
    cromosomaMutado[gen1] = individuo[gen2]
    cromosomaMutado[gen2] = individuo[gen1]
    return cromosomaMutado

def OperadorElitismo(población,puntuaciones):
    '''
    Calcula el índice dentro de la población del cromosoma con mayor puntuación.
    Entradas:
    población: nxm matrix: matriz que contiene los cromosomas
    puntuaciones: 1xn numpy array: lista de puntuaciones de cada cromosoma
    Salidas:
    imax: int: índice dentro de la población del mejor cromosoma
    '''
    imax = 0
    for i in range(len(población)): #se encuentra el mejor cromosoma
        if puntuaciones[i]>puntuaciones[imax]:
            imax = i
    return imax

def GraficarRuta(cromosoma,ciudadesXY,generación,subdirectorio):
    '''
    Crea una gráfica con matplotlib de la trayectoria seguida por la secuencia de un cromosoma.
    Entradas:
    cromosoma: 1xn matrix: cromosoma a graficar
    ciudadesXY: 2xn matrix: matriz [X,Y] que contiene las coordenadas de cada ciudad
    generación: int: ciclo principal de la simulación
    subdirectorio: str: nombre de la carpeta donde se guardará la gráfica
    Salidas:
    Esta función no retorna ningún valor, pero crea un archivo png en el directorio indicado.
    '''
    secuenciasub = np.copy(cromosoma)
    secuencia = np.append(secuenciasub,secuenciasub[0])
    totalGenes = len(secuencia)
    for gen in range(totalGenes-1):
        ciudad1 = secuencia[gen]
        ciudad2 = secuencia[gen+1]
        x1, x2 = ciudadesXY[0][ciudad1], ciudadesXY[0][ciudad2]
        y1, y2 = ciudadesXY[1][ciudad1], ciudadesXY[1][ciudad2]
        xstep = abs(x2-x1)/1000
        ystep = abs(y2-y1)/1000
        if x1<x2:
            x = np.arange(x1,x2,xstep)
            if y1<y2:
                y = np.arange(y1,y2,ystep)
            elif y1>y2:
                y = np.flip(np.arange(y2,y1,ystep))
            else:
                y = np.ones(1000)
                y[:]*=y1
        elif x1>x2:
            x = np.arange(x2,x1,xstep)
            if y1<y2:
                y = np.flip(np.arange(y1,y2,ystep))
            elif y1>y2:
                y = np.arange(y2,y1,ystep)
            else:
                y = np.ones(1000)
                y[:] *= y1
        else:
            x = np.ones(1000)
            x[:] *= x1
            if y1<y2:
                y = np.arange(y1,y2,ystep)
            else:
                y = np.flip(np.arange(y2,y1,ystep))
        plt.plot(x[0:1000],y[0:1000],'b')
    plt.plot(ciudadesXY[0],ciudadesXY[1],'r*')
    genes = len(cromosoma)
    plt.plot(ciudadesXY[0][(cromosoma[0])],ciudadesXY[1][(cromosoma[0])],'g*',label='inicio')
    plt.plot(ciudadesXY[0][(cromosoma[genes-1])],ciudadesXY[1][(cromosoma[genes-1])],'y*',label='fin')
    plt.title('Mejor ruta en gen: '+str(generación))
    endPath = subdirectorio+'/'+'Ruta en gen '+str(generación)+'.png'
    plt.savefig(endPath)
    plt.close('all')
    return

def LeerArchivo(nombreTxt):
    '''
    Obtiene los puntos en XY de las ciudades a utilizar a partir de una matriz almacenada en un archivo
    txt. Los retorna como una matriz [X,Y]
    Entradas:
    nombreTxt: str: nombre del archivo .txt que se va a leer
    Salidas:
    [X,Y]: 2xn matrix: matriz [X,Y] que contiene las coordenadas de cada ciudad
    '''
    datosTxt = open(nombreTxt,'rt')
    listaDeStr = (datosTxt.read()).split('\n')
    for m in range(len(listaDeStr)):
        string = listaDeStr[m]
        listaDeStr[m] = string.replace(']','').replace('[','').replace(' ','')
    for l in range(len(listaDeStr)):
        string = listaDeStr[l]
        listaDeStr[l] = string[0:(len(string)-1)]
    listaDeStr = listaDeStr[0:(len(listaDeStr)-1)]
    X = []
    Y = []
    for x in range(len(listaDeStr)):
        X.append(float(listaDeStr[x].split(',')[0]))
        Y.append(float(listaDeStr[x].split(',')[1]))
    return [X,Y]

def EscribirArchivo(matriz,nombreTxt):
    '''
    Escribe un archivo txt de la matriz que se indica en el mismo directorio que este script.
    Entradas:
    matriz: numpy array: array que se debe almacenar en el archivo txt
    nombreTxt: str: nombre del archivo .txt que se va a almacenar
    Salidas:
    Esta función no retorna ningún valor, pero crea un archivo txt en el directorio del script.
    '''
    with open(nombreTxt,'w') as archivo:
        archivo.write(str(matriz))
    return

def Optimización():
    '''
    Algoritmo principal del AGE. Los parámetros se indican en las zonas apropiadas. La función
    no retorna ningún valor, pero llama a funciones que crean archivos txt y png, junto con sus
    directorios necesarios.
    '''
    np.random.seed(23432)
    # Parámetros iniciales
    tamañoPoblación = 20
    nGeneraciones = 10000
    nombreTxt = 'CoordenadasCiudades.txt'
    ubicaciónCiudades = LeerArchivo(nombreTxt)
    nGenes = len(ubicaciónCiudades[0])
    # Inicialización de variables
    subfolder = CrearDirectorio('Gráficas AGE')
    mejorPuntuación = 0
    mejoresLongitudes = []
    población = InicializarPoblación(tamañoPoblación, nGenes)
    generaciónvector = range(nGeneraciones)
    LongitudesProm = []
    for generación in range(nGeneraciones): #cada ciclo es una generación completa
        puntuaciones = EvaluarPoblación(población,ubicaciónCiudades) #se evalúa a la población
        icromosomaMax = OperadorElitismo(población,puntuaciones) #se selecciona el que no va a cambiar
        mejoresLongitudes.append(1/puntuaciones[icromosomaMax]) #se guarda la mejor longitud actual
        LongitudesProm = np.append(LongitudesProm,1/np.average(puntuaciones)) #se guarda la puntuación promedio
        if puntuaciones[icromosomaMax]>mejorPuntuación: #se grafica la nueva mejor ruta
            mejorPuntuación = puntuaciones[icromosomaMax]
            GraficarRuta(población[icromosomaMax],ubicaciónCiudades,generación,subfolder)
        for i in range(tamañoPoblación): #se modifica a todos menos el mejor
            if i!=icromosomaMax:
                población[i] = OperadorMutación(población[i])
    EscribirArchivo(población[icromosomaMax],'caminoMásCorto_AGE.txt') #se guarda la ruta más corta
    GraficarRuta(población[icromosomaMax],ubicaciónCiudades,nGeneraciones,subfolder) #se grafica la última mejor ruta
    print('Longitud mínima: '+str(mejoresLongitudes[len(mejoresLongitudes)-1]))
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(generaciónvector, mejoresLongitudes, 'b',label = 'Longitud más corta')
    ax.plot(generaciónvector, LongitudesProm, 'r',label = 'Longitud promedio')
    ax.set_title('Longitud vs generación')
    ax.legend()
    fig.savefig('LvsGen AGE')
    plt.close('all')
    return

Optimización()
