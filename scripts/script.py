from pathlib import Path
from playwright.sync_api import sync_playwright, TimeoutError
import time

# Pasta para salvar os downloads
DOWNLOAD_DIR = Path(r"D:\downloads\proteinas")
DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

# URLs e seletores
BROWSE_URL = "https://www.sasbdb.org/browse/"
PROTEIN_LINK_SELECTOR = "a[href^='/data/']"
DOWNLOAD_BUTTON_SELECTOR = "#download_menu_title"
FULL_ENTRY_LINK_SELECTOR = "a:has-text('full entry (zip)')"

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        print(f"Iniciando na página de navegação: {BROWSE_URL}...")
        page.goto(BROWSE_URL, timeout=0)

        # --- FASE 1: Coletar TODOS os links de TODAS as páginas ---
        all_links = []
        page_num = 1
        while True:
            print(f"Coletando links da página {page_num}...")
            
            page.wait_for_selector(PROTEIN_LINK_SELECTOR, timeout=0)
            
            result_elems = page.query_selector_all(PROTEIN_LINK_SELECTOR)
            links_on_page = ["https://www.sasbdb.org" + el.get_attribute("href") for el in result_elems]
            all_links.extend(links_on_page)

            next_button = page.locator('a:has-text("Next")')
            
            if next_button.count() == 0:
                print("Última página alcançada.")
                break
            
            # <<< CORREÇÃO FINAL: LÓGICA DE ESPERA POR MUDANÇA DE TEXTO >>>
            # 1. Antes de clicar, salvamos o TEXTO do primeiro link da tabela ATUAL.
            first_link_text = page.locator(PROTEIN_LINK_SELECTOR).first.inner_text()
            
            print("Indo para a próxima página...")
            # -------------------------------------------------------------
            # AQUI ESTÁ A CORREÇÃO: timeout=0 e no_wait_after=True
            # -------------------------------------------------------------
            next_button.click(timeout=0, no_wait_after=True)
            
            # 2. ESPERA DEFINITIVA: Usamos uma função JS para esperar até que o texto do 
            #    primeiro link da tabela seja DIFERENTE do texto da página antiga.
            print("Aguardando o texto do conteúdo ser atualizado...")
            page.wait_for_function(f"""
                () => document.querySelector("a[href^='/data/']").innerText !== "{first_link_text}"
            """, timeout=0)
            
            page_num += 1

        all_links = list(dict.fromkeys(all_links))
        print(f"\nColeta finalizada. Total de {len(all_links)} links únicos encontrados em {page_num} páginas.")
        print("Iniciando fase de downloads...")
        time.sleep(2)
        # --- FIM DA FASE 1 ---

        # --- FASE 2: Baixar os arquivos da lista completa ---
        for i, link in enumerate(all_links, start=1):
            print(f"\rProcessando [{i}/{len(all_links)}]: {link}", end="", flush=True)
            try:
                page.goto(link, timeout=0)
                download_button = page.locator(DOWNLOAD_BUTTON_SELECTOR)
                if download_button.count() == 0:
                     print(f"\n  -> Seletor '#download_menu_title' não encontrado. Pulando.", flush=True)
                     continue
                
                download_button.click()
                full_entry_link = page.locator(FULL_ENTRY_LINK_SELECTOR)
                full_entry_link.wait_for(state="visible", timeout=10000)

                with page.expect_download() as download_info:
                    full_entry_link.click()
                
                download = download_info.value
                save_path = DOWNLOAD_DIR / download.suggested_filename
                download.save_as(str(save_path))

                print(f"\r{' ' * 120}\r", end="") 
                print(f"✅ [{i}/{len(all_links)}] Baixado: {save_path.name}", flush=True)

            except TimeoutError:
                print(f"\n  -> ⚠️  Timeout no link {i}. O item de download pode não existir. Pulando.", flush=True)
            except Exception as e:
                print(f"\n  -> ❌ Erro inesperado no link {i}: {e}", flush=True)
            finally:
                time.sleep(0.2)

        print("\n\nProcesso finalizado!")
        browser.close()

if __name__ == "__main__":
    run()
