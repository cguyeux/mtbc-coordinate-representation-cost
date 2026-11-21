"""Client scriptable du serveur MCP TBannotator (transport HTTP streamable).

Permet à un script (cron, orchestrateur) d'exécuter une requête SQL en lecture
seule sur la base TBannotator et de récupérer le résultat en DataFrame, sans
dépendre de l'agent Claude. Endpoint et transport repris de la config du harness
(`~/.claude.json` : tbannotator = https://tblearn.tbannotator.ideev.universite-paris-saclay.fr/mcp).

Nécessite le SDK `mcp` (installé dans le venv `.venv-ingest`). Le serveur tronque
au-delà de `max_rows` ; `query_df` pagine automatiquement par LIMIT/OFFSET si besoin.
"""
import io

import anyio
import pandas as pd
from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client

MCP_URL = "https://tblearn.tbannotator.ideev.universite-paris-saclay.fr/mcp"
PAGE = 50000


async def _call(sql: str, max_rows: int, timeout_s: int):
    async with streamablehttp_client(MCP_URL) as (read, write, _):
        async with ClientSession(read, write) as s:
            await s.initialize()
            r = await s.call_tool("tool_query_postgres",
                                  {"query": sql, "max_rows": max_rows,
                                   "timeout_seconds": timeout_s})
            # le serveur renvoie un repr de dict Python ; parse robuste via ast
            txt = r.content[0].text
            obj = _loads_pyrepr(txt)
            if not obj.get("success"):
                raise RuntimeError(obj.get("error", txt[:300]))
            d = obj["data"]
            return d["csv"], d.get("row_count", 0), d.get("truncated", False)


def _loads_pyrepr(txt: str) -> dict:
    """Le MCP renvoie un repr de dict Python (quotes simples). Parse via ast."""
    import ast
    return ast.literal_eval(txt)


def query_df(sql: str, max_rows: int = PAGE, timeout_s: int = 60,
             paginate: bool = True) -> pd.DataFrame:
    """Exécute `sql` et renvoie un DataFrame. Pagine si le résultat est tronqué.

    Pour paginer, `sql` ne doit pas déjà contenir LIMIT/OFFSET ; on les ajoute.
    """
    csv0, n0, trunc = anyio.run(_call, sql, max_rows, timeout_s)
    df = pd.read_csv(io.StringIO(csv0)) if csv0.strip() else pd.DataFrame()
    if not (trunc and paginate):
        return df
    # pagination
    frames, offset = [df], len(df)
    base = sql.rstrip().rstrip(";")
    while True:
        page_sql = f"{base} LIMIT {PAGE} OFFSET {offset}"
        csv_i, n_i, trunc_i = anyio.run(_call, page_sql, PAGE, timeout_s)
        di = pd.read_csv(io.StringIO(csv_i)) if csv_i.strip() else pd.DataFrame()
        if di.empty:
            break
        frames.append(di)
        offset += len(di)
        if len(di) < PAGE:
            break
    return pd.concat(frames, ignore_index=True)
