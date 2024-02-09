import { defineConfig } from 'vite'
import { svelte } from '@sveltejs/vite-plugin-svelte'

// https://vitejs.dev/config/
export default defineConfig({
  bas: "/spdc_cavity_analysis/",
  plugins: [svelte()],
})
