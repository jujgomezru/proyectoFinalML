"""Seleccion reproducible de muestras de PRJEB11755 a partir de una consulta en vivo a ENA.

Se pide a la API de ENA el listado completo de corridas del proyecto (590 en total) y se descarta
la corrida "companera" (chica, <1GB) de cada muestra biologica -- cada una tiene una corrida
profunda (~2-5GB, la que se procesa) y una companion insuficiente para ensamblar -- lo que deja
295 corridas profundas repartidas en los grupos del proyecto (alias por granja/cohorte: PIG, EYZ,
SYZ, BMZ, BHZ, ZXZ, DB).

`seleccionar_runs(n_muestras=None)` devuelve las 295 corridas profundas completas -- el universo
real a procesar. Con `n_muestras` en un numero, en cambio, reparte ese numero entre los grupos
proporcionalmente a `PESOS_GRUPO`, con el metodo de mayores restos y al menos 1 muestra por grupo,
prefiriendo dentro de cada grupo las corridas mas livianas -- un subconjunto reducido para
prototipar el pipeline sin pagar el costo de las 295 (asi lo usa `04_pipeline_organizado.ipynb`).

Solo se consideran corridas `library_layout=PAIRED` -- el pipeline de descarga (`resolver_y_descargar()`
en `02_estimacion_recursos.ipynb`/`04_pipeline_organizado.ipynb`) asume siempre dos archivos FASTQ
(R1/R2) y no reconoce `SINGLE`.

`PESOS_GRUPO` queda fijo (no se recalcula contra el tamano actual de cada grupo en ENA) para que
correr esto de nuevo no reordene ni reemplace las muestras ya elegidas -- y por lo tanto ya
descargadas y procesadas -- en corridas anteriores.
"""
import pandas as pd

PESOS_GRUPO = {"PIG": 200, "EYZ": 52, "SYZ": 38, "BMZ": 34, "BHZ": 24, "ZXZ": 22, "DB": 20}


def _reparto_mayores_restos(n_total, pesos):
    total_peso = sum(pesos.values())
    crudo = {g: n_total * w / total_peso for g, w in pesos.items()}
    asignado = {g: int(v) for g, v in crudo.items()}
    restante = n_total - sum(asignado.values())
    por_resto = sorted(pesos, key=lambda g: crudo[g] - asignado[g], reverse=True)
    for g in por_resto[:restante]:
        asignado[g] += 1
    return asignado


def seleccionar_runs(n_muestras=50, proyecto="PRJEB11755", min_bytes_corrida_profunda=1_000_000_000):
    """Consulta ENA en vivo y devuelve run_accessions de `proyecto`.

    `n_muestras=None` devuelve las 295 corridas profundas completas (todo el universo a procesar,
    sin subselecccion). Con un numero, reparte ese numero entre grupos proporcionalmente a
    `PESOS_GRUPO` (ver modulo) -- pensado para un subconjunto reducido de prototipado.
    """
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={proyecto}&result=read_run"
        "&fields=run_accession,sample_alias,fastq_bytes,library_layout&format=tsv"
    )
    ena = pd.read_csv(url, sep="\t")
    ena["grupo"] = ena["sample_alias"].str.extract(r"^([A-Za-z]+)")
    ena["bytes"] = ena["fastq_bytes"].apply(lambda v: sum(int(x) for x in str(v).split(";")))
    profundas = ena[(ena["bytes"] > min_bytes_corrida_profunda) & ena["grupo"].isin(PESOS_GRUPO)
                     & (ena["library_layout"] == "PAIRED")]

    if n_muestras is None:
        profundas = profundas.sort_values(["grupo", "bytes"])
        return list(profundas["run_accession"])

    asignado = _reparto_mayores_restos(n_muestras, PESOS_GRUPO)
    runs = []
    for grupo, n in asignado.items():
        del_grupo = profundas[profundas["grupo"] == grupo].sort_values("bytes")
        runs.extend(del_grupo["run_accession"].head(n))
    return runs
